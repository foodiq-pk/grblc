import os
from abc import ABC
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from stsci.tools import capable

from config import CONFIG
from external.PythonPhot import aper
from grblc.data_processing.datastructures import Image, SkyObject


class TransformatorManager:
    """transformator manager for applying transformators to whole image list
    Takes transform list in constructor and applies transformation onto passed image list"""
    def __init__(self, transform_list):
        transform_list.sort(key=lambda x: x.VALUE)
        self.transformator = Transformator(transform_list)

    def apply_transformations(self, image_list):
        transformed_list = []
        for image in image_list:
            transformed_list.append(self.transformator.apply(image))
        return transformed_list


class Transformator:
    """
    Class to combine all transforms to be applied on an image.
    """

    def __init__(self, transformators):
        """
        Tranformators to be applied on image.
        Must be initialized in order of application.
        :param transformators: list of transformators
        """
        transformators.sort(key=lambda x: x.VALUE)
        self.transformators = transformators

    def apply(self, image: Image):
        """
        applies all initialized transformators on given image
        :param image: Image
        :return: modified image by all given transformators
        """
        cpimage = deepcopy(image)
        for transformator in self.transformators:
            cpimage = transformator.transform(cpimage)
        return cpimage


class BaseTransform(ABC):

    REQUIRES = []
    PROVIDES = []

    def transform(self, image: Image):
        """
        Applies transformation to an image
        :param image: Image
        :return: modified Image / does not copy
        """
        raise NotImplementedError

    def requirements_check(self, image: Image):
        """
        Checks for necessary prerequisites for given transformation.
        :param image: Image
        :return: False if image is missing any of the required transformation
        """
        return False not in [req in image.processing_parameters for req in self.REQUIRES]


class FlatTransform(BaseTransform):
    """
    Transform that applies flat field correction.

    Must be initialized with master flat Image object or with create_master_flat function to create one
    Gives option to create master flat image from flat images using static method "create_master_flat"

    """
    REQUIRES = ["dark"]
    PROVIDES = ["flat"]
    VALUE = 2

    def __init__(self, master_flat: Image):
        try:
            image = fits.open(master_flat.fixed_parameters["path"])
            self.flat = image[0].data
        except AttributeError:
            raise AttributeError("Accepts only single Image type object as an argument")

    @staticmethod
    def create_master_flat(images: [Image], save_path: str =CONFIG["FLAT_PATH"]):
        """
        Static method to create master flat from given flat images. Save path specified in Config.FLAT_PATH
        Overwrites flat image in Config.FLAT_PATH !
        :param images: Flat images list
        :param save_path: Optional save path
        :return: master flat as Image type object
        """
        flat_filename = "master_flat" + datetime.now().strftime("_%H-%M-%Y-%d-%m") + ".fits"
        save_path = Path(save_path) / flat_filename
        counter = 0
        values = []
        for image in images:
            values.append(fits.open(image.fixed_parameters["path"])[0].data)
            counter += 1
        flat_data = np.median(values, axis=0)
        flat_data = flat_data / np.mean(flat_data)
        hdu = fits.PrimaryHDU(flat_data)
        if os.path.exists(save_path) and CONFIG["OVERWRITE"]:
            os.remove(save_path)
        hdu.writeto(save_path)
        return Image({"time_jd": 0, "exposure": 0, "type": "flat", "path": save_path, "id": "mflat"}, {})

    def transform(self, image: Image):
        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        if "flat" in image.processing_parameters:
            raise AttributeError("flat correction already applied")
        img = fits.open(image.fixed_parameters["path"])
        values = img[0].data/self.flat
        header = img[0].header

        # path and filename
        old_path = image.fixed_parameters["path"]
        new_path = old_path.parent
        new_file_name = old_path.stem + "df.fits"
        new_path = new_path / new_file_name
        image.fixed_parameters["path"] = new_path

        # writing fits and creating image object
        if os.path.exists(new_path) and CONFIG["OVERWRITE"]:
            os.remove(new_path)
        fits.writeto(image.fixed_parameters["path"], values, header)
        image.processing_parameters["flat"] = True
        return image


class DarkTransform(BaseTransform):
    """
    Transform that subtracts dark current from image.

    Must be initialized with master dark Image object or with create_master_dark function to create one
    Gives option to create master flat frame using static method "create_master_dark"
    """
    REQUIRES = []
    PROVIDES = ["dark"]
    # value to order transformations in a logicla order of application
    VALUE = 1

    def __init__(self, master_dark: Image):
        try:
            self.image = master_dark
            img = fits.open(master_dark.fixed_parameters["path"])
            self.dark = img[0].data
        except AttributeError:
            raise AttributeError("Accepts only single Image type object as an argument")

    @staticmethod
    def create_master_dark(images: [Image], exposure=20,  save_path=CONFIG["DARK_PATH"]):
        """
        Static method to create master flat from given flat images. Save path specified in Config.DARK_PATH

        :param images: Dark images list
        :param save_path: Optional save path
        :param exposure: exposure of input dark images to be selected and combined
        :return: master flat as Image type object
        """
        dark_filename = "master_dark" + datetime.now().strftime("_%H-%M-%Y-%d-%m") + ".fits"
        save_path = Path(save_path) / dark_filename
        counter = 0
        values = []
        for image in images:
            if image.get_exposure() == exposure:
                values.append(fits.open(image.fixed_parameters["path"])[0].data)
                counter += 1
        if counter == 0:
            raise RuntimeError("No darks with given exposure")
        dark_data = np.median(values, axis=0)
        hdu = fits.PrimaryHDU(dark_data)
        if os.path.exists(save_path) and CONFIG["OVERWRITE"]:
            os.remove(save_path)
        hdu.writeto(save_path)
        return Image({"time_jd": 0, "exposure": 20, "type": "dark", "path": save_path, "id": "mdark"}, {})

    def transform(self, image: Image):

        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        if "dark" in image.processing_parameters:
            raise AttributeError("dark correction already applied")
        img = fits.open(image.fixed_parameters["path"])
        values = img[0].data - image.get_exposure() / self.image.get_exposure() * self.dark
        header = img[0].header


        # path and filename
        old_path = image.fixed_parameters["path"]
        new_path = old_path.parent
        new_file_name = old_path.stem + "d.fits"
        new_path = new_path / new_file_name
        image.fixed_parameters["path"] = new_path
        # overrwriting old one
        if os.path.exists(new_path) and CONFIG["OVERWRITE"]:
            os.remove(new_path)

        # writing image to fits and creating image object
        fits.writeto(image.fixed_parameters["path"], values, header)
        image.processing_parameters["dark"] = True
        return image


@DeprecationWarning
class PyrafPhotometryTransform(BaseTransform):
    """
    Transform that does daophot tast from IRAF on image with objects specified in inicialization.

    Must be initialized with object list of Object type.
    Writes result in processed parameteres as photometry entry.
    ** TODO:
        - Writing results to objects
    """
    REQUIRES = ["dark", "flat"]
    PROVIDES = ["photometry"]
    # value to order transformations in a logicla order of application
    VALUE = 3

    def __init__(self, object_list: [SkyObject]):
        self.objects = object_list

    def _get_coordinates_file(self, image: Image):
        """ Creates coordinates file for photometry using wcs coordinates in FITS header.

        :return: filename
        """
        filename = CONFIG["FILE_DUMP"] + "xy_coords_file.txt"
        header = fits.getheader(image.get_path())
        w = WCS(header)
        pixel_coordinates = []

        for o in self.objects:
            x, y = w.wcs_world2pix(o.fixed_parameters['ra'],
                                   o.fixed_parameters['dec'], 1)

            pixel_coordinates.append([float(x),
                                      float(y)])
        np.savetxt(filename, pixel_coordinates)
        return filename

    def transform(self, image: Image):

        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        # INIT
        from pyraf import iraf
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.daophot(_doprint=0)
        capable.OF_GRAPHICS = False

        # File handling
        temp_final_out = CONFIG["FILE_DUMP"] + 'TempPhotOut.dat'
        if os.path.exists(temp_final_out):
            os.remove(temp_final_out)

        base = image.fixed_parameters["path"].split('/')[-1]
        dao_phot_out = CONFIG["FILE_DUMP"] + base + ".mag.dat"
        phot_txt_out = CONFIG["FILE_DUMP"] + base + "PhoTxOut.dat"

        # Setting pyraf phot parameters
        photpars = iraf.photpars.getParList()
        iraf.photpars.setParam('apertures', CONFIG["APERTURE"])
        iraf.phot.setParam('image', image.fixed_parameters["path"])
        iraf.phot.setParam('coords', self._get_coordinates_file(image))
        iraf.phot.setParam('verify', 'no')
        iraf.phot.setParam('output', dao_phot_out)
        iraf.phot.setParam('interactive', 'no')

        if os.path.exists(dao_phot_out):
            os.remove(dao_phot_out)

        # Running IRAF task
        dump = iraf.phot(mode='h', Stdout=1)

        if os.path.exists(phot_txt_out):
            os.remove(phot_txt_out)

        # Getting better formatted file with magnitudes
        iraf.txdump(dao_phot_out, 'XCENTER,YCENTER,MAG,MERR', 'yes', Stdout=phot_txt_out)

        # Getting results from output file
        with open(phot_txt_out) as output_file:
            results = {}
            i = 0
            for lines in output_file:
                parts = lines.split()
                results[self.objects[i].fixed_parameters["id"]] = (float(parts[2]) if parts[2] != "INDEF" else None,
                                                                   float(parts[3]) if parts[3] != "INDEF" else None)
                i += 1

        image.processing_parameters["photometry"] = results
        return image


class PythonPhotPhotometryTransform(BaseTransform):
    """
    PythonPhot photometry procedure

    Must be initialized with object list of Object type.
    Writes result in processed parameteres as photometry entry.
    ** TODO:
        - Writing results to objects
    """
    REQUIRES = ["dark", "flat"]
    PROVIDES = ["photometry"]
    # value to order transformations in a logicla order of application
    VALUE = 3

    def __init__(self, object_list: [SkyObject], aperture=CONFIG['APERTURE']):
        self.objects = object_list
        self.aperture = aperture

    def _get_coordinates_arrays(self, image: Image):
        """ Creates coordinates arrays for photometry using wcs coordinates in FITS header.

        :return: filename
        """
        header = fits.getheader(image.get_path())
        w = WCS(header)
        pixel_coordinates_ra = []
        pixel_coordinates_dec = []

        for o in self.objects:
            x, y = w.wcs_world2pix(o.fixed_parameters['ra'],
                                   o.fixed_parameters['dec'], 1)

            pixel_coordinates_ra.append(float(x))
            pixel_coordinates_dec.append(float(y))
        return np.array(pixel_coordinates_ra), np.array(pixel_coordinates_dec)

    def transform(self, image: Image):
        img = fits.getdata(image.get_path())
        xpos, ypos = self._get_coordinates_arrays(image)
        mag, magerr, flux, fluxerr, sky, skyerr, badflag, outstr = \
            aper.aper(img, xpos, ypos, phpadu=1, apr=self.aperture, zeropoint=25,
                      skyrad=[40, 50], badpix=[-12000, 1060000],
                      exact=True)
        results = {}
        i = 0
        for o in self.objects:

            results[o.get_id()] = (float(mag[i]) if mag[i] is not "nan" else None,
                                   float(magerr[i]) if magerr[i] is not "nan" else None)
            i += 1

        image.processing_parameters["photometry"] = results
        # workaround for single object photometry
        try:
            image.processing_parameters["src_flux"] = (float(flux[0]), float(fluxerr[0]))
            image.processing_parameters["sky"] = (float(sky[0]), float(skyerr[0]))
        except IndexError:
            image.processing_parameters["src_flux"] = (float(flux), float(fluxerr))
            image.processing_parameters["sky"] = (float(sky), float(skyerr))
        return image


class ShiftTransform(BaseTransform):
    """
    Transform to calculate difference in zero magnitude of a frame.
    Must be initialized with object list that is the same as object list used for photometry on given image to calculate
    shifts of each of the control stars.

    Writes information about shift and also specific values for each of the stars in the processing parameters under
    shift and shifts respectively.
    """
    REQUIRES = ["photometry"]
    PROVIDES = ["shifts"]

    # value to order transformations in a logicla order of application
    VALUE = 4

    def __init__(self, object_list: [SkyObject]):
        ol = []
        for sobject in object_list:  # removing grbs from list to calc shifts
            if sobject.get_type() == "star":
                ol.append(sobject)
        self.object_list = ol

    def find_star_with_id(self, id):
        for star in self.object_list:
            if id == star.get_id():
                return star

        return None

    @staticmethod
    def calculate_shift(image):
        vals = image.get_shifts().values()
        return ShiftTransform.weighted_mean(vals)

    @staticmethod
    def weighted_mean(pairs):
        sig_pow = 2
        sigma2 = 1 / (np.sum([1 / s[1] ** sig_pow for s in pairs]))
        mean = np.sum([s[0] / s[1] ** sig_pow for s in pairs]) * sigma2
        return mean, np.sqrt(sigma2)

    def transform(self, image: Image):
        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        shifts = {}
        for id, value in image.get_photometry().items():
            star = self.find_star_with_id(id)
            if star is None:
                continue
            cat_mag = star.get_catalog_magnitude()
            shifts[id] = (cat_mag[0]-value[0],
                          np.sqrt(cat_mag[1]**2+value[1]**2))

        image.processing_parameters["shifts"] = shifts
        image.processing_parameters["shift"] = ShiftTransform.calculate_shift(image)
        return image
