from abc import ABC
import pyfits
from data_tier_placeholder.image import Image
import numpy as np


class BaseTransform(ABC):

    def transform(self, image: Image):
        raise NotImplementedError


class FlatTransform(BaseTransform):
    """
    Transform that applies flat image correction.

    Requires image to be corrected by dark images.

    Gives option to create master flat image from flat images using static method "create_master_flat"

    """
    REQUIRES = ["dark"]
    PROVIDES = ["flat"]

    def __init__(self, master_flat: Image):
        image = pyfits.open(master_flat.fixed_parameters["path"])
        self.flat = image[0].data

    @staticmethod
    def create_master_flat(images: [Image], save_path: str ="/tmp/temp_master_flat.fits"):
        """
        Static method to create master flat from given flat images, with option to save output, default save is /tmp/
        :param images: Flat images list
        :param save_path: Optional save path
        :return: master flat as Image type object
        """
        counter = 0
        values = []
        for image in images:
            values.append(pyfits.open(image.fixed_parameters["path"]))
            counter += 1
        flat_data = np.median(values, axis=0)
        flat_data = flat_data / np.mean(flat_data)
        hdu = pyfits.PrimaryHDU(flat_data)
        hdu.writeto(save_path)
        return Image({"time_jd": 0, "exposure": 0, "type": "flat", "path": save_path})

    def transform(self, image: Image):
        pass
        # TODO


class DarkTransform(BaseTransform):
    """
    Transform that applies dark frame correction.
    Gives option to create master flat frame using static method "create_master_dark"
    """
    REQUIRES = []
    PROVIDES = ["dark"]

    def __init__(self):
        pass
        # TODO

    @staticmethod
    def create_master_dark(images: [Image], exposure: float):
        pass
        # TODO

    def transform(self, image: Image):
        pass
        # TODO


class PhotometryTransform(BaseTransform):

    def __init__(self):
        pass
        # TODO

    def transform(self, image: Image):
        pass
        # TODO


class ShiftTransform(BaseTransform):

    def __init__(self):
        pass
        # TODO

    def transform(self, image: Image):
        pass
        # TODO


class StackTransform(BaseTransform):

    def __init__(self):
        pass
    # TODO

    def transform(self, image: Image):
        pass
        # TODO


class CalibrateTransform(BaseTransform):

    def __init__(self):
        pass
        # TODO

    def transform(self, image: Image):
        pass
        # TODO
