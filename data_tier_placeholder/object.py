from data_tier_placeholder.datastructureabc import DataStructure


class Object(DataStructure):
    """

    """

    REQUIRES = ["ra", "dec", "catalog_magnitude", "id"]
    # TODO difference between GRB and stars

    def __init__(self, fixed_parameters: dict, processing_parameters: dict = {}):
        super().__init__(fixed_parameters, processing_parameters)
