from matplotlib import pyplot as plt

from grblc.data_processing.datahandlers import *
from grblc.image_processing.transformators import *
from grblc.image_processing.transformators import Transformator


def plot_src_flux(image_list: [Image]):
    time = []
    src_flux = []

    for img in image_list:
        time.append(img.get_time_jd())
        src_flux.append(img.get_src_flux()[0])
    plt.plot(time, src_flux, "x")


def plot_sn(image_list: [Image]):
    time = []
    sn = []

    for img in image_list:
        time.append(img.get_time_jd())
        sn.append(img.get_src_flux()[0]/ img.get_src_flux()[1])
    plt.plot(time,sn,"+")


def plot_sky(image_list: [Image]):
    time = []
    sky = []
    for img in image_list:
        time.append(img.get_time_jd())
        sky.append(img.get_sky()[0])
    plt.plot(time, sky, "x")


if __name__ == "__main__":
    # TODO: make grb accesible for plotting later with seconds from trigger as it is more valuable information
    # IS MAINLY TO SHOW HOW TO USE SCRIPT, NOT YET BUILT FOR SERIOUS USE.

    # TODO: recreate config for use and changing parameters there.
    # set up some configuration beforehand
    # Config.set_aperture(4)
    # Config.set_dark_path('/tmp')
    # Config.set_flat_path('/tmp')
    # Config.set_db("sqlite:///tmp/temp.db")
    # OVERWRITE = True
    # DATA_DIR = "/home/foodiq/data/grbs/"
    # # dialog window
    # Config.print()

    # get all frames using dialog windows
    # darks = FileHandlerDialog("dark").get_list()
    # flats = FileHandlerDialog("flat").get_list()
    # data = FileHandlerDialog("data").get_list()

    # get all frames using FileHandler
    darks = FileHandler(folder="/home/foodiq/data/grbs/131030A/d50/c0/darks/",
                        query="*RA.fits",
                        data_type="dark").get_list()[:10]
    flats = FileHandler(folder="/home/foodiq/data/grbs/131030A/d50/c0/flats/",
                        query="*RA.fits",
                        data_type="flat").get_list()[:10]

    data = FileHandler(folder="/home/foodiq/data/grbs/131030A/d50/c0/51018/",
                       query="*RA.fits",
                       data_type="data").get_list()[:10]

    # initialize grb object
    grb = SkyObject({"ra": 345.06729,
                     "dec": -5.3684,
                     "jd_trigger": 2456596.372431,
                     "type": "grb", "id": "GRB131030A",
                     "catalog_magnitude": (None, None)}, {})

    # get objects in vicinity for calibration
    object_list = ObjectHandler(grb).get_list()

    # create master dark from darks loaded
    mdark = DarkTransform.create_master_dark(darks)

    # correct flat images (create transformator manager with dark transform and apply to all flats
    dark_transformator_man = Transformator([DarkTransform(mdark)])

    # TODO: change to one method (image_list, transformator)
    corrected_flats = []
    for flat_image in flats:
        corrected_flats.append(dark_transformator_man.apply(flat_image))

    # create master flat
    mflat = FlatTransform.create_master_flat(corrected_flats)

    # create all transformators - dark and flat correction, photometry and shift
    photometry_transform_manager = Transformator([DarkTransform(mdark),
                                                 FlatTransform(mflat),
                                                 PythonPhotPhotometryTransform(object_list),
                                                 ShiftTransform(object_list)])

    # apply to all data images
    output_data = []
    for image in data:
        output_data.append(photometry_transform_manager.apply(image))

    # save results to database file specified in config
    # overwriting just now, gonna change later

    # plot data
    object_list[0].plot_light_curve(output_data)

    # load data from database (loads from db specified in config or any other path specified when creating
    #  DatabaseHandler)
    # loaded_images, loaded_objects = DatabaseHandler(Config.DB_ENGINE).get_list()

    # loaded_objects[0].plot_lightcurve(loaded_images)

