from data_tier_placeholder.datastructureabc import DataStructure


class Object(DataStructure):
    """
    # TODO
    """

    REQUIRES = ["ra", "dec", "catalog_magnitude" or "trigger_jd", "id", "type"]
    # TODO difference between GRB and stars

    def __init__(self, fixed_parameters: dict, processing_parameters: dict = {}):
        super().__init__(fixed_parameters, processing_parameters)
