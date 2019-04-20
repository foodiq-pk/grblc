import os
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import ccdproc
from astropy.io import fits
from astropy.wcs import WCS
from ccdproc import wcs_project

from config import CONFIG
from grblc.data_processing.datastructures import Image
from grblc.image_processing.transformators import PythonPhotPhotometryTransform, Transformator


class StackingManager:
    """class for handling stacking procedures"""
    pass
    # TODO


def stack_until_err_or_limit(image_list, limiting_magnitude_error, grb, max_frames=15):
    """
    Stacks images when it runs into image where magnitude error is over limit, tries to stack with next frame,
    does photometry again, if it passes through limitin error goes to next frame and repeats or adds another image.
    Maximum amount of images in one stack (through addition of one per cycle of stack + photometry + check)
    can be specified
    :param image_list: list of Image type objects sorted by time of exposure ascending
    :param limiting_magnitude_error: error for passing image as okay without need to stack more or at all
    :param max_frames: limitin amount of frames so you get better coverage for lightcurve
    :return: image list with stacked images instead of original ones. Sorted ascending by time
    """
    # create temp folder for this stacking procedure with timestamp (for later cleanup)
    time_init = datetime.now().strftime("%Y_%m_%d_T%H-%M-%S")
    folder_path = Path("/tmp/stack" + time_init)
    # TODO: check  or separate to fction???
    folder_path.mkdir()

    # filter input
    image_list = _simple_filter_bad_images(image_list, s_n=0.5)

    phot = Transformator([PythonPhotPhotometryTransform([grb])])
    output_image_list = []
    i = 0
    print("len = ", len(image_list))
    while i < len(image_list):
        print("i: ", i)
        image = image_list[i]
        s_n = image.get_src_flux()[0]/image.get_src_flux()[1]
        # test this one
        names = []
        if "stack" in image.processing_parameters:
            names = image.get_stack
        else:
            pass
        # TODO: check the condition (S/N or error)

        if s_n < limiting_magnitude_error:
            # get some image info for stacked image
            time_start = image.get_time_jd()

            # for time separation limitation in the future
            # TODO:
            exposure = image.get_exposure()

            stacked_image = image
            header = fits.getheader(stacked_image.get_path())
            wcs_header = WCS(header)
            ccddata_repro_stack = ccdproc.CCDData(data=fits.getdata(stacked_image.get_path()),
                                                  wcs=wcs_header
                                                  , unit="adu")
            reprojected = [ccddata_repro_stack]

            n = 1
            print("pred checkem sn")
            while s_n < limiting_magnitude_error or n > max_frames:
                print("sn = ", s_n)
                n += 1

                # TODO: add time separation condition
                if n+i >= len(image_list):
                    print("ende")
                    # no more images to stack or images time separation is over limit
                    output_image_list.append(stacked_image)
                    i += n
                    break
                else:
                    # get image to add info
                    print("i+n", i+n)
                    image_to_add = image_list[i+n]
                    data = fits.getdata(image_to_add.get_path())
                    header_wcs = fits.getheader(image_to_add.get_path())

                    # create ccd data type image of the one to add
                    ccddata = ccdproc.CCDData(data, wcs=WCS(header_wcs), unit="adu")

                    # add the image after reprojecting onto first image WCS
                    reprojected.append(wcs_project(ccddata, wcs_header))

                    # stack images
                    combiner = ccdproc.Combiner(reprojected)
                    stacked_image_data = combiner.average_combine()
                    stacked_image_data.wcs = WCS(header)

                    # get last image data
                    time_end = image_to_add.get_time_jd() + image_to_add.get_exposure() / 86000
                    time_jd_mid = (time_start + time_end) / 2
                    time_coverage = (time_end - time_start) * 86000
                    # modify header and create new filename
                    header_wcs["EXPOSURE"] = time_coverage
                    header["JD"] = time_jd_mid
                    filename = "stack-" + str(time_jd_mid) + "-e" + "{:.0f}".format(time_coverage) + ".fits"

                    # TODO: saving of stacks for later use? garbage cleaning? (finished list difference with real
                    #  images in a folder - delete ones that are not in output list after stacking)
                    # temporary save path

                    file_path = folder_path / filename

                    # new file shenanigans
                    if os.path.exists(file_path) and CONFIG["OVERWRITE"]:
                        os.remove(file_path)
                    fits.writeto(file_path, stacked_image_data, header)
                    names.append(image_to_add.get_path())

                    # new image object
                    stacked_image = Image(fixed_parameters={"path": file_path,
                                                            "exposure": time_coverage,
                                                            "time_jd": time_jd_mid,
                                                            "type": "data", "id": filename},
                                          processing_parameters={"flat": True,
                                                                 "dark": True,
                                                                 "stack": names})
                    # do phot on grb and get s/n of stack to check for improvement
                    stacked_image = phot.apply(stacked_image)
                    s_n = stacked_image.get_src_flux()[0] / stacked_image.get_src_flux()[1]
            # TODO: check for duplicates (last image twice??)
            print("saved: " + str(stacked_image.get_path()))
            output_image_list.append(stacked_image)
            i +=n + 1

        else:
            # if magnitude or S/N ratio is sufficient append to output without stacking
            output_image_list.append(image)
            i += 1

    return output_image_list


def stack_two_neighbours(image_list, limiting_magnitude_error, grb, max_reruns=2):
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
    # TODO: test
    def pick_images_for_combining_and_pop_them(image_list: [Image], ratio_limit):
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
        # print(len(images_to_stack), len(image_list))
        return images_to_stack, image_list

    def stack_and_append_to_old(images_to_stack: [[Image]], image_list: [Image]):
        # print("img flagger for stacking in colors:")
        # plot_stacks(grb, images_to_stack=images_to_stack, image_list=image_list)
        # print("stacking")
        for images in images_to_stack:
            image_list.append(stacking_procedure(images))
        return image_list

    def phot_redo_temp(image_list: [Image]):
        # phot trans only for grb
        phot = Transformator([PythonPhotPhotometryTransform(grb)])

        image_list_phot = []
        # print("phot")
        for img in image_list:
            # print("phot for {}".format(img.get_path()))
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


def stack_using_sequence(image_list, sequence):
    # TODO: TEST
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
    final_list.sort(key = (lambda x: x.get_time_jd()))
    return final_list


def _simple_filter_bad_images(image_list, s_n):
    """
    throw away images that are below rough limiting magnitude precision. Can be e.g. 1mag
    error or so depending on time of observation after burst trigger.
    :param image_list:
    :param limiting_magnitude_error:
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
            names.append(image.get_path())
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
    header["JD"] = timejd_mid

    filename = "stack-" + str(timejd_mid) + "-e" + "{:.0f}".format(timecoverage) + ".fits"

    path = Path("/tmp")
    path = path / filename

    if os.path.exists(path):
        os.remove(path)
    fits.writeto(path, final_image_data, header)

    stacked_image = Image(fixed_parameters={"path": path,
                                            "exposure": timecoverage,
                                            "time_jd": timejd_mid,
                                            "type": "data", "id": filename},
                          processing_parameters={"flat": True,
                                                 "dark": True,
                                                 "stack": names})
    return stacked_image
