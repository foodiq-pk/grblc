from abc import ABC
import pyfits
from data_tier_placeholder.image import Image
import numpy as np
import os


class BaseTransform(ABC):

    REQUIRES = []
    PROVIDES = []

    def transform(self, image: Image):
        raise NotImplementedError

    def requirements_check(self, image: Image):
        return False not in [req in image.processing_parameters for req in self.REQUIRES]


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
        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        img = pyfits.open(image.fixed_parameters["path"])
        values = img[0].data/self.flat
        header = img[0].header
        old_path = img.fixed_parameters["path"]
        img.fixed_parameters["path"] = old_path[:5] + "f.fits"
        pyfits.writeto(img.fixed_parameters["path"], values, header)
        os.remove(old_path)
        image.processing_parameters["flat"] = True


class DarkTransform(BaseTransform):
    """
    Transform that applies dark frame correction.
    Gives option to create master flat frame using static method "create_master_dark"
    """
    REQUIRES = []
    PROVIDES = ["dark"]

    def __init__(self, master_dark: Image):
        image = pyfits.open(master_dark.fixed_parameters["path"])
        self.dark = image[0].data

    @staticmethod
    def create_master_dark(images: [Image], save_path='/tmp/temp_master_dark.fits'):
        """
        Static method to create master dark from given flat images, with option to save output, default save is /tmp/
        :param images: Dark images list
        :param save_path: Optional save path
        :return: master flat as Image type object
        """
        counter = 0
        values = []
        for image in images:
            values.append(pyfits.open(image.fixed_parameters["path"]))
            counter += 1
        dark_data = np.median(values, axis=0)
        hdu = pyfits.PrimaryHDU(dark_data)
        hdu.writeto(save_path)
        return Image({"time_jd": 0, "exposure": 0, "type": "dark", "path": save_path})

    def transform(self, image: Image):
        if not self.requirements_check(image):
            raise ValueError("Missing required transformations on image. Need:" + str(self.REQUIRES))
        img = pyfits.open(image.fixed_parameters["path"])
        values = img[0].data - self.dark
        header = img[0].header
        old_path = img.fixed_parameters["path"]
        img.fixed_parameters["path"] = old_path[:5] + "d.fits"
        pyfits.writeto(img.fixed_parameters["path"], values, header)
        image.processing_parameters["dark"] = True


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
