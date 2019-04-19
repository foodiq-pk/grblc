from data_tier_placeholder.datastructureabc import DataStructure
from astropy.io import fits
from astropy.wcs import WCS


class Image(DataStructure):
    """
    Image class data object, contains
    """
    REQUIRES = ["time_jd", "exposure", "path", "type", "id"]

    def __init__(self, fixed_parameters: dict, processing_parameters: dict = {}):
        super().__init__(fixed_parameters, processing_parameters)

    def print_info(self):
        print(self.fixed_parameters,
              self.processing_parameters)

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
