from data_tier_placeholder.image import Image
from astropy.io import fits
import ccdproc
from astropy.wcs import WCS
from ccdproc import wcs_project
import os
from pathlib import Path


class StackingManager:
    """class for handling stacking procedures"""
    pass
    # TODO


def _stack_until_err_or_limit(image_list, limiting_magnitude_error, max_frames):
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
    pass
    # TODO


def _stack_two_neighbours(image_list, limiting_magnitude_error, max_reruns):
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
    pass
    # TODO


def _stack_using_specified_sequence(image_list, sequence):
    """
    Stacks images according to specified sequence, sequence is a list of numbers saying how many frames it should stack
    next. Example:
        [1,1,1,1,1,1,2,2,5,5,5]
        This sequence would not stack first six images then stack 7th and 8th, 9th and 10th, and 11th to 16th,
         17th to 22nd etc.
    Sequence should be long enought to cover lenth of image list.
    If there is a frame that is too far timewise from other frames (Time mid difference of frames
    is more than twice their exposure e.g.) stack frames until that point and continue with frames after the gap with
    the same number from sequence and continue after the required amount is stacked at least once.
    :param image_list: list of Image type objects sorted ascending by time of exposure
    :param sequence: list of numbers specifying amounts of frames stacked in each cycle while going through image list
    :return: image list with stacked images sorted ascending by time
    """
    pass
    # TODO


def _simple_filter_bad_images(image_list, limiting_magnitude_error):
    """
    throw away images that are below rough limiting magnitude precision. Can be e.g. 1mag
    error or so depending on time of observation after burst trigger.
    :param image_list:
    :param limiting_magnitude_error:
    :return:
    """
    return [img for img in image_list if img.get_photometry()[0][1] < limiting_magnitude_error]


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
                                            "type": "cdata", "id": filename},
                          processing_parameters={"flat": True,
                                                 "dark": True,
                                                 "stack": names})
    return stacked_image
