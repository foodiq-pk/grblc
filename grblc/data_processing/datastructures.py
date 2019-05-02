from abc import ABC

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt


class DataStructure(ABC):
    """
    Abstract class for data structures
    """
    REQUIRES = []

    def __init__(self, fixed_parameters: dict, processing_parameters: dict = {}):
        """
        constructor for both data structures
        :param fixed_parameters: permanent characteristics
        :param processing_parameters: Parameters calculated from images with transformators
        """

        if self.check_required_input_parameters(fixed_parameters):
            self.fixed_parameters = fixed_parameters
        else:
            raise ValueError("Missing one or more required parameters.")
        self.processing_parameters = processing_parameters

    def check_required_input_parameters(self, parameters: dict):
        """
        Checks for legal input in image constructor.
        :param parameters: parameter dictionary
        :return: bool, whether input is legal
        """
        return all([req in parameters for req in self.REQUIRES])

    def get_type(self):
        return self.fixed_parameters["type"]

    def get_id(self):
        return self.fixed_parameters["id"]


class Image(DataStructure):
    """
    Image class data object, contains
    """
    REQUIRES = ["time_jd", "exposure", "path", "type", "id"]

    def __init__(self, fixed_parameters: dict, processing_parameters: dict = {}):
        super().__init__(fixed_parameters, processing_parameters)
        # TODO: logging debug creation of this (repr?)

    def __repr__(self):
        return f"Image('{self.get_id()}'," \
               f"'{self.get_time_jd():16.8f}'," \
               f"'{self.get_exposure()}'," \
               f"'{self.get_type()}'," \
               f"'{self.processing_parameters.keys()}')"

    def get_time_jd(self):
        return self.fixed_parameters["time_jd"]

    def get_exposure(self):
        return self.fixed_parameters["exposure"]

    def get_path(self):
        return self.fixed_parameters["path"]

    def get_photometry(self):
        return self.processing_parameters["photometry"]

    def get_shifts(self):
        return self.processing_parameters["shifts"]

    def get_shift(self):
        return self.processing_parameters["shift"]

    def get_src_flux(self):
        return self.processing_parameters["src_flux"]

    def get_sky(self):
        return self.processing_parameters["sky"]

    def get_stack(self):
        return self.processing_parameters["stack"]

    def check_if_object_are_in_image(self, object_list):
        """
        Checks whether objects in passed objects list are within image borders or not.

        :return: tuple, list of objects in image and list of objects not in image
        """
        header = fits.getheader(self.get_path())
        data = fits.getdata(self.get_path())
        shape = data.shape
        x_lim = shape[1]
        y_lim = shape[0]
        print(shape)
        w = WCS(header)
        in_image = []
        not_in_image = []

        for o in object_list:
            x, y = w.wcs_world2pix(o.fixed_parameters['ra'],
                                   o.fixed_parameters['dec'], 1)
            if x < x_lim and y < y_lim:
                in_image.append(o)
            else:
                not_in_image.append(o)
        return in_image, not_in_image


class SkyObject(DataStructure):
    """
    # TODO
    """

    REQUIRES = ["ra", "dec", "id", "type"]

    def __init__(self, fixed_parameters: dict, processing_parameters: dict = {}):
        super().__init__(fixed_parameters, processing_parameters)
        # TODO: logging debug creation of this (repr?)

    def __repr__(self):
        return f"SkyObject('{self.get_id()}'," \
               f"'{self.get_ra()}'," \
               f"'{self.get_dec()}'," \
               f"'{self.get_type()}'," \
               f"'{self.processing_parameters.keys()}')"

    @staticmethod
    def star(id: str, ra: float, dec: float, catalog_magnitude: tuple, cat_filter: str):
        """
        Static method to create a star object
        :param id:
        :param ra:
        :param dec:
        :param catalog_magnitude: tuple with catalog magnitudes
        :return:
        """
        # TODO: docs
        fixed_pars = {"id": id,
                      "ra": ra,
                      "dec": dec,
                      "type": "star"}
        processing_pars = {"catalog_magnitude": catalog_magnitude,
                           "catalog_filter": cat_filter}

        return SkyObject(fixed_parameters=fixed_pars, processing_parameters=processing_pars)

    @staticmethod
    def grb(name, ra, dec, trigger_jd):
        """
        Static method to create grb object
        :param name:
        :param ra:
        :param dec:
        :param trigger_jd:
        :return:
        """
        # TODO: docs

        fixed_pars = {"id": name,
                      "ra": ra,
                      "dec": dec,
                      "type": "grb"
                      }
        processing_pars = {"jd_trigger": trigger_jd}
        return SkyObject(fixed_parameters=fixed_pars, processing_parameters=processing_pars)

    def get_ra(self):
        return self.fixed_parameters["ra"]

    def get_dec(self):
        return self.fixed_parameters["dec"]

    def get_catalog_magnitude(self):
        if self.get_type() == "star":
            return self.processing_parameters["catalog_magnitude"]
        else:
            return None, None

    def get_catalog_filter(self):
        if self.get_type() == "star":
            return self.processing_parameters["catalog_filter"]
        else:
            return None, None

    def get_trigger_jd(self):
        if self.get_type() == "grb":
            return self.processing_parameters["jd_trigger"]
        else:
            raise TypeError("Star object doesnt have JD Trigger")

    def get_raw_light_curve(self, image_list):
        times = []
        time_errs = []
        values = []
        value_errs = []

        for img in image_list:
            mag = img.get_photometry()[self.get_id()][0]
            magerr = img.get_photometry()[self.get_id()][1]
            if mag is not None:
                times.append(img.get_time_jd() + img.get_exposure()/84600/2)
                time_errs.append(img.get_exposure()/84600/2)
                values.append(mag)
                value_errs.append(magerr)
        return times, time_errs, values, value_errs

    def get_shifted_light_curve(self, image_list):
        times = []
        time_errs = []
        values = []
        value_errs = []

        for img in image_list:
            mag = img.get_photometry()[self.get_id()][0]
            magerr = img.get_photometry()[self.get_id()][1]
            if mag is not None:
                times.append(img.get_time_jd() + img.get_exposure()/84600/2)
                time_errs.append(img.get_exposure()/84600/2)
                values.append(mag + img.get_shift()[0])
                value_errs.append( np.sqrt(magerr**2 + img.get_shift()[1]**2))
        return times, time_errs, values, value_errs

    def plot_light_curve(self, image_list, type="raw", timeerr=False, magerr=False):
        if type is "raw":
            time, time_errs, mag, mag_errs = self.get_raw_light_curve(image_list)
        elif type is "shifted":
            time, time_errs, mag, mag_errs = self.get_shifted_light_curve(image_list)
        else:
            raise ValueError("Invalid light curve type")

        if timeerr is False:
            time_errs = None
        if magerr is False:
            mag_errs = None
        plt.errorbar(time, mag, yerr=mag_errs, xerr=time_errs, fmt=".")
        plt.xlabel("Time - JD[days]")
        plt.ylabel("Magnitude [mag]")
        plt.grid()
        plt.gca().invert_yaxis()
        plt.show()
