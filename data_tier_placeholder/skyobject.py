from data_tier_placeholder.datastructureabc import DataStructure
import numpy as np


class SkyObject(DataStructure):
    """
    # TODO
    """

    REQUIRES = ["ra", "dec", "id", "type"]  # TODO add restriction "catalog_magnitude" or "trigger_jd"

    def __init__(self, fixed_parameters: dict, processing_parameters: dict = {}):
        super().__init__(fixed_parameters, processing_parameters)

    def get_ra(self):
        return self.fixed_parameters["ra"]

    def get_dec(self):
        return self.fixed_parameters["dec"]

    def get_catalog_magnitude(self):
        return self.fixed_parameters["catalog_magnitude"]

    def get_raw_light_curve(self, image_list):
        times = []
        values = []
        value_errs = []

        for img in image_list:
            mag = img.get_photometry()[self.get_id()][0]
            magerr = img.get_photometry()[self.get_id()][1]
            if mag is not None:
                times.append(img.get_time_jd())
                values.append(mag)
                value_errs.append(magerr)
        return times, values, value_errs

    def get_shifted_light_curve(self, image_list):
        times = []
        values = []
        value_errs = []

        for img in image_list:
            mag = img.get_photometry()[self.get_id()][0]
            magerr = img.get_photometry()[self.get_id()][1]
            if mag is not None:
                times.append(img.get_time_jd())
                values.append(mag - img.get_shift()[0])
                value_errs.append( np.sqrt(magerr**2 + img.get_shift()[1]**2))
        return times, values, value_errs
