from data_tier_placeholder.datastructureabc import DataStructure


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

