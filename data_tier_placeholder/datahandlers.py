import glob
from abc import ABC
from pathlib import Path
from tkinter import Tk
from tkinter.filedialog import askopenfilenames

import astropy.coordinates as coord
import astropy.units as u
from astropy.io import fits as pyfits
from astroquery.vizier import Vizier
from sqlalchemy.engine import create_engine
from sqlalchemy.orm import sessionmaker

from data_tier_placeholder.database import Frame, Magnitude, Shift, SObject, init_session, Base
from data_tier_placeholder.image import Image
from data_tier_placeholder.skyobject import SkyObject
from config import Config

# TODO: check for possible crashes because of this change (fits to pyfits)


class InvalidQueryError(Exception):
    pass


class DataManagerFactory:
    """
    Factory class for selecting desired handlers based on path and query type.
    """
    @staticmethod
    def get_handler(path: str, query: str=None):
        """
        Supplies handler based on type of path:
        "sqlite:" gets database handler with path to it given db file (e.g.  sqlite:/tmp/temp.db)
        "/tmp/"   gets file handler with getting all files in a folder restricted by query
        "/tmp/file.fits" gets file handler for single file
        query format: shell regex with space delimiter for additional type specification
        :param path:
        :param query: a query string in format "regex type_specification" for example "*df.fits data"
        :return: a handler for type of path and query
        """
        pass
        # TODO add allowed type options (data/flat/dark/  possibly - grb/star)


class BasicHandler(ABC):

    def __init__(self, path: str, query: str=None):
        self.query = query
        self.path = Path(path)

    def get_list(self):
        """returns list of objects based on selected handler and query"""
        raise NotImplementedError


class FileHandler(BasicHandler):
    # TODO: add docs
    """
    Handler for file and folder type queries.
    Only allows generation of lists of images and single image objects
    Examples of usage:

    """
    def __init__(self, folder: str, query=None, data_type="data"):
        super().__init__(folder, query)
        if data_type not in ["flat", "dark", "data", "cdata", "dfdata"]:
            raise Exception("Invalid data type specified")

        self.type = data_type
        if self.query is None:
            self.selector = "*.fits"
        else:
            self.selector = query  # TODO: CHECK FOR VALIDITY OF SELECTOR

    def get_list(self):
        # data type and corrections
        if self.type in ["flat", "dark", "data"]:
            processing_pars = {}
        elif self.type in ["cdata", "dfdata"]:
            processing_pars = {"flat": True, "dark": True}
            self.type = "data"
        else:
            processing_pars = {"dark": True}
            self.type = "data"

        p = self.path
        paths = list(p.glob(self.selector))
        image_list = []
        for im_path in paths:
                image = pyfits.open(im_path)
                exposure = image[0].header["EXPOSURE"]
                time_jd = image[0].header["JD"]
                fixed_pars = {
                    "time_jd": time_jd,
                    "path": im_path,
                    "type": self.type,
                    "exposure": exposure,
                    "id": im_path.stem}
                image_list.append(Image(fixed_parameters=fixed_pars,
                                        processing_parameters=processing_pars))

        image_list.sort(key=lambda x: x.get_time_jd())
        return image_list


class FileHandlerDialog:
    """DIalog version of file handler for interactive usage"""
    def __init__(self, type_string):
        if type_string in ["flat", "dark", "data", "ddata", "cdata", "dfdata"]:
            self.type = type_string
        else:
            raise Exception("Invalid FITS file type")

        self.pathlist = self._gen_path_list_from_dialog()

    def _gen_path_list_from_dialog(self):
        Tk().withdraw()
        files = askopenfilenames(title="choose {0} type FITS files".format(self.type),
                                 initialdir=Config.DATA_DIR,
                                 filetypes=[("FITS files", "*.fits")],
                                 )
        paths = []
        for file in files:
            paths.append(Path(file))

        return paths

    def get_list(self):
        if self.type in ["cdata", "dfdata"]:
            processing_parameters = {"flat": True, "dark": True}
        elif self.type is "ddata":
            processing_parameters = {"dark": True}
        else:
            processing_parameters = {}
        image_list = []
        for path in self.pathlist:
            image = pyfits.open(path)
            exposure = image[0].header["EXPOSURE"]
            time_jd = image[0].header["JD"]
            fixed_pars = {
                "time_jd": time_jd,
                "path": path,
                "type": self.type,
                "exposure": exposure,
                "id": path.stem}

            image_list.append(Image(fixed_parameters=fixed_pars,
                                    processing_parameters=processing_parameters))

        image_list.sort(key=lambda x: x.get_time_jd())
        return image_list


class DatabaseHandler(BasicHandler):
    """
    Handler for interacting with database (sqlite3)
    Allows loading and saving images and object data from database file.
    """
# TODO: redo using pathlib
    def __init__(self, path: str, query: str=None):
        super().__init__(path, query)

    def get_list(self):
        # TODO: WORKAROUND engine path sqlite:///
        engine = create_engine("sqlite:///" + str(self.path))
        _session = sessionmaker()
        session = _session(bind=engine)
        image_list = self._load_images(session)
        object_list = self._load_objects(session)
        self._load_photometry(session, image_list)
        self._load_shifts(session, image_list)
        return image_list, object_list

    @staticmethod
    def _load_images(session):
        image_list = []
        db_image_list = session.query(Frame).all()
        for db_image in db_image_list:
            image_list.append(Image(fixed_parameters={"time_jd": db_image.time_jd,
                                                      "exposure": db_image.exposure,
                                                      "id": db_image.id,
                                                      "type": db_image.type,
                                                      "path": Path(db_image.path)},
                                    processing_parameters={"dark": True,
                                                           "flat": True,
                                                           "shift": (db_image.shift, db_image.shifterr)}))
        return image_list

    @staticmethod
    def _load_objects(session):

        object_list = []
        db_obj_list = session.query(SObject).all()
        for db_obj in db_obj_list:
            object_list.append(SkyObject(fixed_parameters={"ra": db_obj.ra,
                                                           "dec": db_obj.dec,
                                                           "id": db_obj.id,
                                                           "type": "star",
                                                           "catalog_magnitude": (db_obj.catmag, db_obj.catmagerr)}))
        return object_list

    @staticmethod
    def _load_photometry(session, image_list):
        for image in image_list:
            star_list = session.query(Magnitude).filter(Magnitude.frame_id == image.fixed_parameters["id"]).all()
            photometry = {}
            for star in star_list:
                photometry[star.star_id] = (star.mag, star.magerr)
            image.processing_parameters["photometry"] = photometry

    @staticmethod
    def _load_shifts(session, image_list):
        for image in image_list:
            star_list = session.query(Shift).filter(Shift.frame_id == image.fixed_parameters["id"]).all()
            shifts = {}
            for star in star_list:
                shifts[star.star_id] = (star.shift, star.shifterr)
            image.processing_parameters["shifts"] = shifts

    @staticmethod
    def save_objects_and_images(image_list: [Image], object_list: [SkyObject]):
        """Saves everything from object list and image list to database specified in Config.DB_ENGINE.
        Destination db can be changed through session engine parameter or in config"""

        session = init_session()

        # frames
        for image in image_list:
            session.add(Frame(id=image.get_id(),
                              time_jd=image.get_time_jd(),
                              type=image.get_type(),
                              exposure=image.get_exposure(),
                              path=str(image.get_path()),
                              shift=image.get_shift()[0],
                              shifterr=image.get_shift()[1]))

        # objects
        for skyobject in object_list:
            session.add(SObject(id=str(skyobject.get_id()),
                                ra=skyobject.get_ra(),
                                dec=skyobject.get_dec(),
                                catmag=skyobject.get_catalog_magnitude()[0],
                                catmagerr=skyobject.get_catalog_magnitude()[1]))

        # processing pars
        for image in image_list:
            for star_id, magnitude in image.get_photometry().items():
                session.add(Magnitude(star_id=str(star_id),
                                      frame_id=image.get_id(),
                                      mag=magnitude[0],
                                      magerr=magnitude[1]))
            # TODO: formatting of Shifts parameter after calculating (for saving-missing ID-not needed only for tests)
        for image in image_list:
            for star_id, shift in image.get_shifts().items():
                session.add(Shift(star_id=str(star_id),
                                  frame_id=image.get_id(),
                                  shift=shift[0],
                                  shifterr=shift[1]))
        Base.metadata.create_all()
        session.commit()
    # TODO: location of base, engine etc.


class ObjectHandler:
    """Handler to create list of objects in specified vicinity of given target

    Supported catalogues:
        APASS
        NOMAD

    Usage:
        If not specified the default values for area around target coordinates are:
            radius = 0.1 degrees
            catalog: APASS
            lower magnitude limit: 16 mag
            maximum objects is 100

        Initialize with target object (GRB) and specify result limit if necessary. You can change default values in
        get_list() call as follows:

            ObjectHandler(target).get_list(mag_limit = 20, catalog= "NOMAD", radius = 2 * u.deg)

         * - Units for radius are in degrees

    """
    # TODO Static methods?

    def __init__(self, target: SkyObject, result_limit=100):
        self.target = target
        self.limit = result_limit

    def get_list(self, mag_limit=16., catalog="APASS", radius=0.1,):
        """

        :param mag_limit: float
        :param catalog: catalog name string
        :param radius: radius in degrees
        :return: Object type list with target and objects in specified vicinity
        """
        if catalog == "APASS":
            return self.vizier_query_object_list_apass(self.target, radius, mag_limit)
        if catalog == "NOMAD":
            return self.vizier_query_object_list_nomad(self.target, radius, mag_limit)
        else:
            raise ValueError("Unsupported or invalid catalog name")

    def vizier_query_object_list_apass(self, target: SkyObject, radius, mag_limit):
        """
        Creates object list from vizier APASS query

        :param target: Center point around which to query for object
        :param radius: radius of area around central object in degrees
        :param mag_limit: limiting lower magnitude for star selection
        :return: list of objects and the GRB
        """
        vizier_query = Vizier(columns=['recno', 'RAJ2000', 'DEJ2000', 'Vmag', 'e_Vmag'],
                              column_filters={"Vmag": "<"+str(mag_limit)},
                              row_limit=self.limit)
        coordinates = coord.SkyCoord(target.fixed_parameters["ra"],
                                     target.fixed_parameters["dec"],
                                     unit=(u.deg, u.deg),
                                     frame='icrs')
        result = vizier_query.query_region(coordinates,
                                           radius=radius * u.deg,
                                           catalog="APASS")

        object_list = [target]
        i = 1
        for o in result[0]:
            object_list.append(SkyObject({"ra": o['RAJ2000'],
                                          "dec": o['DEJ2000'],
                                          "catalog_magnitude": (o['Vmag'], o['e_Vmag']),
                                          "id": o["recno"],
                                          "type": "star"}))
            i += 1

        return object_list

    def vizier_query_object_list_nomad(self, target: SkyObject, radius, mag_limit):
        """
        Creates object list from vizier NOMAD query

        :param target: Center point around which to query for object
        :param radius: radius of area around central object in degrees
        :param mag_limit: limiting lower magnitude for star selection
        :return: list of objects and the GRB
        """
        vizier_query = Vizier(columns=['RAJ2000', 'DEJ2000', 'Vmag'],
                              column_filters={"Vmag": "<"+str(mag_limit)},
                              row_limit=self.limit)
        coordinates = coord.SkyCoord(target.fixed_parameters["ra"],
                                     target.fixed_parameters["dec"],
                                     unit=(u.deg, u.deg),
                                     frame='icrs')
        result = vizier_query.query_region(coordinates,
                                           radius=radius * u.deg,
                                           catalog="NOMAD")

        object_list = [target]
        i = 1
        for o in result[0]:
            object_list.append(SkyObject({"ra": o['RAJ2000'],
                                          "dec": o['DEJ2000'],
                                          "catalog_magnitude": (o['Vmag'], o['e_Vmag']),
                                          "id": i,
                                          "type": "star"}))
            i += 1

        return object_list
