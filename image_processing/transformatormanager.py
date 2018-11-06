from data_tier_placeholder.image import Image


class Transformator:
    """


    """

    def __init__(self, transformators):
        self.transformators = transformators

    def apply(self, image: Image):
        for transformator in self.transformators:
            image = transformator.transform(image)
        return image

