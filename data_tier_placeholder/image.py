from data_tier_placeholder.datastructureabc import DataStructure


class Image(DataStructure):
    """
    Image class data object, contains
    """
    REQUIRES = ["time_jd", "exposure", "path", "type"]

    def __init__(self, fixed_parameters: dict, processing_parameters: dict = {}):
        super().__init__(fixed_parameters, processing_parameters)


