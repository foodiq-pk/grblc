from data_tier_placeholder.datastructureabc import DataStructure
import numpy as np
import matplotlib.pyplot as plt


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
        time_errs = []
        values = []
        value_errs = []

        for img in image_list:
            mag = img.get_photometry()[self.get_id()][0]
            magerr = img.get_photometry()[self.get_id()][1]
            if mag is not None:
                times.append(img.get_time_jd())
                time_errs.append(img.get_exposure()/84600)
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
                times.append(img.get_time_jd())
                time_errs.append(img.get_exposure()/84600)
                values.append(mag - img.get_shift()[0])
                value_errs.append( np.sqrt(magerr**2 + img.get_shift()[1]**2))
        return times, time_errs, values, value_errs

    def plot_light_curve(self, image_list, type="raw",timeerr = False, magerr = False):
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
        plt.errorbar(time,mag,yerr=mag_errs,xerr=time_errs,fmt=".")
        plt.xlabel("Time - JD[days]")
        plt.ylabel("Magnitude [mag]")
        plt.grid()
        plt.gca().invert_yaxis()
        plt.show()
