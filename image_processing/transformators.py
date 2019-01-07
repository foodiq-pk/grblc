import os
from abc import ABC

import numpy as np
import pyfits
from astropy.io import fits
from astropy.wcs import WCS
from stsci.tools import capable

from config import Config
from data_tier_placeholder.image import Image
from data_tier_placeholder.skyobject import SkyObject


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

    def __init__(self, master_flat: Image):
        image = pyfits.open(master_flat.fixed_parameters["path"])
        self.flat = image[0].data

    @staticmethod
    def create_master_flat(images: [Image], save_path: str =Config.FLAT_PATH):
        """
        Static method to create master flat from given flat images. Save path specified in Config.FLAT_PATH
        Overwrites flat image in Config.FLAT_PATH !
        :param images: Flat images list
        :param save_path: Optional save path
        :return: master flat as Image type object
        """
        counter = 0
        values = []
        for image in images:
            values.append(pyfits.open(image.fixed_parameters["path"])[0].data)
            counter += 1
        flat_data = np.median(values, axis=0)
        flat_data = flat_data / np.mean(flat_data)
        hdu = pyfits.PrimaryHDU(flat_data)
        if os.path.exists(Config.FLAT_PATH):
            os.remove(Config.FLAT_PATH)
        hdu.writeto(save_path)
        return Image({"time_jd": 0, "exposure": 0, "type": "flat", "path": save_path, "id": "mflat"}, {})

    def transform(self, image: Image):
        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        if "flat" in image.processing_parameters:
            raise AttributeError("flat correction already applied")
        img = pyfits.open(image.fixed_parameters["path"])
        values = img[0].data/self.flat
        header = img[0].header
        old_path = image.fixed_parameters["path"]
        image.fixed_parameters["path"] = old_path[:-5] + "f.fits"
        pyfits.writeto(image.fixed_parameters["path"], values, header)
        # os.remove(old_path) TODO: delete old?
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

    def __init__(self, master_dark: Image):
        image = pyfits.open(master_dark.fixed_parameters["path"])
        self.dark = image[0].data

    @staticmethod
    def create_master_dark(images: [Image], save_path=Config.DARK_PATH):
        """
        Static method to create master flat from given flat images. Save path specified in Config.DARK_PATH

        :param images: Dark images list
        :param save_path: Optional save path
        :return: master flat as Image type object
        """
        counter = 0
        values = []
        for image in images:
            values.append(pyfits.open(image.fixed_parameters["path"])[0].data)
            counter += 1
        dark_data = np.median(values, axis=0)
        hdu = pyfits.PrimaryHDU(dark_data)
        if os.path.exists(Config.DARK_PATH):
            os.remove(Config.DARK_PATH)
        hdu.writeto(save_path)
        return Image({"time_jd": 0, "exposure": 0, "type": "dark", "path": save_path, "id": "mdark"}, {})

    def transform(self, image: Image):

        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        if "dark" in image.processing_parameters:
            raise AttributeError("dark correction already applied")
        img = pyfits.open(image.fixed_parameters["path"])
        values = img[0].data - self.dark
        header = img[0].header
        old_path = image.fixed_parameters["path"]
        image.fixed_parameters["path"] = old_path[:-5] + "d.fits"
        pyfits.writeto(image.fixed_parameters["path"], values, header)
        image.processing_parameters["dark"] = True
        return image


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

    def __init__(self, object_list: [SkyObject]):
        self.objects = object_list

    def get_coordinates_file(self, image: Image):
        """ Creates coordinates file for photometry using wcs coordinates in FITS header.

        :return: filename
        """
        filename = Config.FILE_DUMP + "xy_coords_file.txt"
        header = fits.getheader(image.fixed_parameters["path"])
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
        temp_final_out = Config.FILE_DUMP + 'TempPhotOut.dat'
        if os.path.exists(temp_final_out):
            os.remove(temp_final_out)

        base = image.fixed_parameters["path"].split('/')[-1]
        dao_phot_out = Config.FILE_DUMP + base + ".mag.dat"
        phot_txt_out = Config.FILE_DUMP + base + "PhoTxOut.dat"

        # Setting pyraf phot parameters
        photpars = iraf.photpars.getParList()
        iraf.photpars.setParam('apertures', Config.APERTURE)
        iraf.phot.setParam('image', image.fixed_parameters["path"])
        iraf.phot.setParam('coords', self.get_coordinates_file(image))
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


class ShiftTransform(BaseTransform):
    """
    Transform to calculate difference in zero magnitude of a frame.
    Must be initialized with object list that is the same as object list used for photometry on given image to calculate
    shifts of each of the control stars.

    Writes information about shift and also specific values for each of the stars in the processing parameters under
    shift and shifts respectively.
    """
    # TODO separate GRB to not be included in the shift calculation
    REQUIRES = ["photometry"]
    PROVIDES = ["shifts"]

    def __init__(self, object_list: [SkyObject]):
        self.object_list = object_list

    def find_star_with_id(self, id):
        for star in self.object_list:
            if id == star.fixed_parameters["id"]:
                return star

        raise ValueError("no star with given id")

    def transform(self, image: Image):
        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        shifts = {}
        # TODO: crash on trying to shift grb - no cat magnitude/or missing in list
        for id, value in image.processing_parameters["photometry"].items():
            cat_mag = self.find_star_with_id(id).fixed_parameters["catalog_magnitude"]
            shifts[id] = (cat_mag[0]-value[0],
                          np.sqrt(cat_mag[1]**2+value[1]**2))

        image.processing_parameters["shifts"] = shifts
        #image.processing_parameters["shift"] = (np.mean([val[0] for val in shifts]),
        #                                       np.sqrt(sum(val[1]**2 for val in shifts)) /
        #                                        len([val[1] for val in shifts]))
        return image


class StackTransform(BaseTransform):

    REQUIRES = ["dark", "flat"]
    PROVIDES = ["stack"]

    def __init__(self):
        pass
    # TODO

    def transform(self, image: Image):
        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        pass
        # TODO


class CalibrateTransform(BaseTransform):

    REQUIRES = ["shift"]
    PROVIDES = ["lightcurve"]  # << ??

    def __init__(self):
        pass
        # TODO

    def transform(self, image: Image):
        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        pass
        # TODO
