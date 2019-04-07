from data_tier_placeholder.image import Image
from astropy.io import fits
import ccdproc
from astropy.wcs import WCS
from ccdproc import wcs_project
import os
from pathlib import Path

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
    timejd_mid = (tfirst + tlast)/ 2
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
