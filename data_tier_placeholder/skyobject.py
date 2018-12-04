from data_tier_placeholder.datastructureabc import DataStructure


class SkyObject(DataStructure):
    """
    # TODO
    """

    REQUIRES = ["ra", "dec", "id", "type"]  # TODO add restriction "catalog_magnitude" or "trigger_jd"

    def __init__(self, fixed_parameters: dict, processing_parameters: dict = {}):
        super().__init__(fixed_parameters, processing_parameters)
