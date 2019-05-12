import os
from copy import deepcopy
from pathlib import Path

import ccdproc
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from ccdproc import wcs_project

from grblc.data_processing.datastructures import Image
from grblc.image_processing.transformators import PythonPhotPhotometryTransform, Transformator


class StackingManager:
    """class for handling stacking procedure"""

    def __init__(self, image_list, grb=None):
        self.images = image_list
        self.grb = grb
        self.to_stack = None
        self.single_images = None
        self.sn_prediction = None
        self.stacked_images = None

    def plot_stack_prediction(self):
        """lets you visualize what images were selected for stacking,
        needs a grb object when creating stacking manager"""

        try:
            if self.grb is None:
                raise RuntimeError("To plot prediction you need to specify grb object when creating stack manager.")
            gtime, gtimerr, gmag, gmagerr = self.grb.get_raw_light_curve(self.images)
            plt.errorbar(x=(np.array(gtime) - self.grb.get_trigger_jd()) * 86400,
                         xerr=np.array(gtimerr) * 86400,
                         y=gmag,
                         yerr=gmagerr,
                         fmt='k.')
            for image_list in self.to_stack:
                gtime, gtimerr, gmag, gmagerr = self.grb.get_raw_light_curve(image_list)
                plt.errorbar(x=(np.array(gtime) - self.grb.get_trigger_jd()) * 86400,
                             xerr=np.array(gtimerr) * 86400,
                             y=gmag,
                             yerr=gmagerr,
                             fmt='.')
            plt.gca().invert_yaxis()
            plt.grid()
            plt.xlabel("time since trigger [s]")
            plt.ylabel("magnitude [mag]")
            # plt.savefig('/home/foodiq/Documents/diplomka/text_prace/images/stack_example.pdf')
            # plt.xscale("log")
            plt.show()
        except AttributeError:
            raise RuntimeError("No images to stack selected. Run select_images_to_stack() first.")


    # def plot_sn_prediction(self):
    #     pass

    def select_images_to_stack(self, sn_limit):
        """"Picks images based on signal to noise specified limit,
         predicts s/n of images to be stacked to pass the limit"""
        i = 0
        # TODO: logging INFO stacking amount

        self.single_images = []
        self.to_stack = []
        self.sn_prediction = []
        while i < len(self.images):
            image = self.images[i]
            if image.get_src_flux()[0] / image.get_src_flux()[1] < sn_limit:
                stack = []
                src_fluxes = [image.get_src_flux()[0], image.get_src_flux()[1] ** 2]
                sn_sq = src_fluxes[0] ** 2 / src_fluxes[1]

                stack.append(image)
                i += 1
                while sn_sq < sn_limit ** 2 and i < len(self.images):
                    image = self.images[i]
                    stack.append(image)
                    src_fluxes[0] += image.get_src_flux()[0]
                    src_fluxes[1] += image.get_src_flux()[1] ** 2
                    sn_sq = src_fluxes[0] ** 2 / src_fluxes[1]

                    i += 1
                if len(stack) > 1:
                    self.to_stack.append(stack)
                    self.sn_prediction.append(np.sqrt(sn_sq))
                else:
                    self.single_images.append(stack[0])
            else:
                self.single_images.append(image)
                i += 1
        # TODO: logging debug what was created and stacks..

    def stack_images(self):
        "Run stacking procedure on prepared stacks object, needs to have selection done first"
        # TODO: logging info statistics
        try:
            self.stacked_images = []
            for ims in self.to_stack:
                self.stacked_images.append(stacking_procedure(ims))
        except TypeError:
            raise RuntimeError("No images to stack selected. Run select_images_to_stack() first.")

    def stack_images_multicore(self, cores=1):
        """allows user to paralelize stacking procedure to be run on specified amount of cores,
        defaults to one which is equivalent to reular stack images method"""
        try:
            from multiprocessing import Pool
            p = Pool(processes=cores)
            res = p.map(stacking_procedure, self.to_stack)
            self.stacked_images = res
        except TypeError:
            raise RuntimeError("No images to stack selected. Run select_images_to_stack() first.")

    def save_stacks(self, folder):
        """ when specifying folder moves results of stacking from temporary folder to the specified folder
        modifies path of Image object in stacked_images attribute to be affiliated with new location"""
        try:
            import shutil
            destination_folder = Path(folder)
            if not (destination_folder.exists() and destination_folder.is_dir()):
                destination_folder.mkdir(parents=True, exist_ok=True)
            new_path_images = deepcopy(self.stacked_images)
            for image in new_path_images:
                old_path = image.get_path()

                new_path = destination_folder / image.get_path().name

                shutil.copy(old_path, str(destination_folder))
                image.fixed_parameters['path'] = new_path
            self.stacked_images = new_path_images
        except TypeError:
            raise RuntimeError("No images stacked to copy. Need to create stack first.")

    def get_sn_predictions(self):
        return self.sn_prediction

    def get_list(self):
        """returns combined and time sorted list of images after stacking"""
        try:
            merged_list = self.stacked_images + self.single_images
            merged_list.sort(key=lambda x: x.get_time_jd())
            return merged_list
        except TypeError:
            raise RuntimeError("Images not stacked. Run stack_images() first.")


def _stack_two_neighbours(image_list, limiting_magnitude_error, grb, max_reruns=2):
    """
    Goes  through image list and if it finds image that has higher magnitude error than limit adds the next one to it
    if it is withing limiting time window. Always stacks only two images next to each other. Can be run multiple times
    to reduce amount of points and increase precision.
    :param image_list: list of Image type objects sorted ascending by time of exposure
    :param limiting_magnitude_error: error of magnitude to initiate stacking of neighbouring frames
    :param max_reruns: amount of times the double stacking procedure will be done if all frames before it dont reach
    limiting magnitude precision
    :return: image list with stacked images sorted ascending by time
    """
    def pick_images_for_combining_and_pop_them(image_list: [Image], ratio_limit):
        """ selects pairs of images and returns image lists to stack and single images"""
        images_to_stack = []
        image_list = deepcopy(image_list)
        image_list.sort(key=lambda x: x.get_time_jd())
        i = len(image_list)
        while i > 0:
            img1 = image_list[-i]
            img2 = image_list[1 - i]
            if img1.get_src_flux()[0] / img1.get_src_flux()[1] < ratio_limit:
                if abs(img1.get_time_jd() - img2.get_time_jd()) * 86000 < (img1.get_exposure() + img2.get_exposure()):
                    images_to_stack.append([image_list.pop(-i), image_list.pop(1 - i)])
                    i -= 2
                    if i == 1:
                        images_to_stack[-1].append(image_list.pop(-1))
                        i -= 1
                        print("baf")
                else:
                    i -= 1
            else:
                i -= 1
        return images_to_stack, image_list

    def stack_and_append_to_old(images_to_stack: [[Image]], image_list: [Image]):
        """recombining the two lists into one final"""

        for images in images_to_stack:
            image_list.append(stacking_procedure(images))
        return image_list

    def phot_redo_temp(image_list: [Image]):
        """redo photometry on a newly stacked images for a rerun"""
        # phot trans only for grb
        phot = Transformator([PythonPhotPhotometryTransform(grb)])

        image_list_phot = []
        for img in image_list:
            image_list_phot.append(phot.apply(img))
        return image_list_phot

    # filter input
    image_list = _simple_filter_bad_images(image_list, s_n=0.05)

    # main body
    for i in range(max_reruns):
        image_list_to_stack, image_list = pick_images_for_combining_and_pop_them(image_list, limiting_magnitude_error)
        image_list = stack_and_append_to_old(images_to_stack=image_list_to_stack, image_list=image_list)
        phot_redo_temp(image_list)

    return image_list


def _stack_using_sequence(image_list, sequence):
    """
    Stacks images according to specified sequence, sequence is a list of numbers saying how many frames it should stack
    next. Example:
        [1,1,1,1,1,1,2,2,5,5,5]
        This sequence would not stack first six images then stack 7th and 8th, 9th and 10th, and 11th to 16th,
         17th to 22nd etc.
    Sequence should be long enough to cover length of image list.
    If there is a frame that is too far time-wise from other frames (Time mid difference of frames
    is more than twice their exposure e.g.) stack frames until that point and continue with frames after the gap with
    the same number from sequence and continue after the required amount is stacked at least once.
    :param image_list: list of Image type objects sorted ascending by time of exposure
    :param sequence: list of numbers specifying amounts of frames stacked in each cycle while going through image list
    :return: image list with stacked images sorted ascending by time
    """

    def group_by_sequence(image_list, sequence):
        single_images = []
        groups = []
        i = 0
        k = 0
        while i < len(image_list) and k < len(sequence):
            if sequence[k] == 1:
                single_images.append(image_list[i])
                k += 1
                i += 1
                print("plop")
            else:
                image_count = sequence[k]
                groups.append(image_list[i:i + image_count])
                k += 1
                i += image_count
        return single_images, groups

    def stack_and_append_to_old(images_to_stack: [[Image]], image_list: [Image]):
        # print("img flagger for stacking in colors:")
        # plot_stacks(grb, images_to_stack=images_to_stack, image_list=image_list)
        # print("stacking")
        for images in images_to_stack:
            image_list.append(stacking_procedure(images))
        return image_list

    single_images, images_to_stack = group_by_sequence(image_list, sequence)
    final_list = stack_and_append_to_old(images_to_stack, image_list)
    final_list.sort(key=(lambda x: x.get_time_jd()))
    return final_list


def _simple_filter_bad_images(image_list, s_n):
    """
    throw away images that are below rough limiting sn.
    :param image_list:
    :param s_n: limiting s/n for bad images
    :return:
    """
    return [img for img in image_list if img.get_src_flux()[0]/img.get_src_flux()[1] > s_n]


def stacking_procedure(images: [Image]):
    """stacks images in list into single one and writes file
    header will be from middle image in the list.
    Time will be between start of first and end of last exposure with exposure length of this difference.

    :param images list of Image list objects to stack
    """

    # go through images and collect info
    names = []
    i = 0
    half = len(images)//2

    tfirst = images[0].get_time_jd()
    tlast = images[-1].get_time_jd() + images[-1].get_exposure()/86000
    timejd_mid = (tfirst + tlast) / 2
    timecoverage = (tlast - tfirst) * 86000  # s

    for image in images:   # get stack names
        if "stack" in image.processing_parameters:
            for name in image.get_stack():
                names.append(name)
        else:
            names.append(str(image.get_path()))
        i += 1
        if i == half:
            header = fits.getheader(image.get_path())

    midpoint_WCS = WCS(header)

    reprojected = []
    # reproject all images onto the middle one
    for image in images:
        data = fits.getdata(image.get_path())
        header_wcs = fits.getheader(image.get_path())
        ccddata = ccdproc.CCDData(data, wcs=WCS(header_wcs), unit="adu")
        reprojected.append(wcs_project(ccddata, midpoint_WCS))

    combiner = ccdproc.Combiner(reprojected)
    final_image_data = combiner.average_combine()
    final_image_data.wcs = WCS(header)
    header["EXPTIME"] = timecoverage
    header["EXPOSURE"] = timecoverage
    header["JD"] = tfirst

    filename = "stack-" + str(timejd_mid) + "-e" + "{:.0f}".format(timecoverage) + ".fits"

    path = Path("/tmp")
    path = path / filename

    if os.path.exists(path):
        os.remove(path)
    fits.writeto(path, final_image_data, header)

    stacked_image = Image(fixed_parameters={"path": path,
                                            "exposure": timecoverage,
                                            "time_jd": tfirst,
                                            "type": "data", "id": filename},
                          processing_parameters={"flat": True,
                                                 "dark": True,
                                                 "stack": names})
    return stacked_image
