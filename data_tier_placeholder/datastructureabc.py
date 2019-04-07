from abc import ABC


class DataStructure(ABC):
    """

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

