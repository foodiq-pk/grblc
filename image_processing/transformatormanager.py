from copy import deepcopy

from data_tier_placeholder.image import Image


class Transformator:
    """
    Class to combine all transforms to be applied on an image.
    """

    def __init__(self, transformators):
        """
        Tranformators to be applied on image.
        Must be initialized in order of application.
        :param transformators: list of transformators
        """
        self.transformators = transformators

    def apply(self, image: Image):
        """
        applies all initialized transformators on given image
        :param image: Image
        :return: modified image by all given transformators
        """
        cpimage = deepcopy(image)
        for transformator in self.transformators:
            cpimage = transformator.transform(cpimage)
        return cpimage

