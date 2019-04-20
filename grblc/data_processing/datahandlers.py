from abc import ABC
from pathlib import Path
from tkinter import Tk
from tkinter.filedialog import askopenfilenames

import astropy.coordinates as coord
import astropy.units as u
from astropy.io import fits
from astroquery.vizier import Vizier
from sqlalchemy.engine import create_engine
from sqlalchemy.orm import sessionmaker

import grblc.data_processing.database.models as db
from config import CONFIG
from grblc.data_processing.datastructures import Image, SkyObject


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
            raise AttributeError("Invalid FITS file type")

        self.type = data_type
        if self.query is None:
            self.selector = "*.fits"
        else:
            self.selector = query  # TODO: CHECK FOR VALIDITY OF SELECTOR
        # TODO: logging debug creation of this (repr?)

    def __repr__(self):
        return f"FileHandler('{self.type}','{str(self.path)}','{self.selector}')"

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
                image = fits.open(im_path)
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
        # TODO: logging debug - images
        # TODO: logging info - len images
        return image_list


class FileHandlerDialog:
    """DIalog version of file handler for interactive usage"""
    def __init__(self, type_string):
        if type_string in ["flat", "dark", "data", "ddata", "cdata", "dfdata"]:
            self.type = type_string
        else:
            raise AttributeError("Invalid FITS file type")

        self.pathlist = self._gen_path_list_from_dialog()
        # TODO: logging debug creation of this (repr?)

    def __repr__(self):
        return f"FileHandlerDialog('{self.type}')"

    def _gen_path_list_from_dialog(self):
        Tk().withdraw()
        files = askopenfilenames(title="choose {0} type FITS files".format(self.type),
                                 initialdir=CONFIG["DATA_DIR"],
                                 filetypes=[("FITS files", "*.fits")],
                                 )
        paths = []
        for file in files:
            paths.append(Path(file))
        # TODO: logging debug - paths prints
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
            image = fits.open(path)
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
        # TODO: logging debug - images
        # TODO: logging info - len images
        return image_list


class DatabaseHandler:
    """
    Handler for interacting with database (sqlite3)
    Allows loading and saving images and object data from database file.
    """
    def __init__(self, connector):
        self.connector = connector
        engine = create_engine(connector)
        db.Base.metadata.bind = engine
        db.Base.metadata.create_all()
        Session = sessionmaker()
        self.session = Session(bind=engine)
        # TODO: logging debug creation of this (repr?)

    def __repr__(self):
        return f"DatabaseHandler('{self.connector}"

    def get_list(self):
        image_list = self._load_images()
        object_list = self._load_objects()
        return image_list, object_list

    def _load_images(self):
        image_list = []
        db_image_list = self.session.query(db.Frame).all()
        # TODO: logging debug - list query
        # TODO: logging info - len of loaded images
        for db_image in db_image_list:
            image_list.append(Image(fixed_parameters={"time_jd": db_image.time_jd,
                                                      "exposure": db_image.exposure,
                                                      "id": db_image.id,
                                                      "type": db_image.type,
                                                      "path": Path(db_image.path)},
                                    processing_parameters=eval(db_image.additional)))
            print(db_image.additional)
        return image_list

    def _load_objects(self):
        object_list = []
        db_obj_list = self.session.query(db.SObject).all()
        # TODO: logging debug - list query
        # TODO: logging info - len of loaded objects
        for db_obj in db_obj_list:
            object_list.append(SkyObject(fixed_parameters={"ra": db_obj.ra,
                                                           "dec": db_obj.dec,
                                                           "id": db_obj.id,
                                                           "type": db_obj.type},
                                         processing_parameters=eval(db_obj.additional)))
        return object_list

    def save_objects_and_images(self, image_list: [Image], object_list: [SkyObject]):
        """

        :param image_list:
        :param object_list:
        :return:
        """
        # frames
        for image in image_list:
            additional = image.processing_parameters
            frame = db.Frame(id=image.get_id(),
                             time_jd=image.get_time_jd(),
                             type=image.get_type(),
                             exposure=image.get_exposure(),
                             path=str(image.get_path()),
                             additional=str(additional))

            self.session.merge(frame)
            # TODO: logging debug added frame
            # TODO: logging info added frames - len(image_list)

        # objects
        for skyobject in object_list:
            additional = skyobject.processing_parameters
            sobject = db.SObject(id=str(skyobject.get_id()),
                                 ra=skyobject.get_ra(),
                                 dec=skyobject.get_dec(),
                                 additional=str(additional))

            self.session.merge(sobject)
            # TODO: logging debug added object
            # TODO: logging info added objects number - len(object list)
        self.session.commit()


class ObjectHandler:
    # TODO update docs
    """Handler to create list of objects in specified vicinity of given target


    """

    def __init__(self, target: SkyObject, result_limit=100):
        self.target = target
        self.limit = result_limit

    def __repr__(self):
        return f"ObjectHandler('{self.target}','{self.limit}')"

    # TODO: logging debug found objects and query pars
    # TODO: logging info amount of objects found in radius cat and mag
    def get_list(self, mag_limit=16., catalog="APASS", radius=0.1,):
        """

        :param mag_limit: float
        :param catalog: catalog name string
        :param radius: radius in degrees
        :return: Object type list with target and objects in specified vicinity
        """
        # TODO: logging debug calling and pars
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
            object_list.append(SkyObject.star(ra=o['RAJ2000'],
                                              dec=o['DEJ2000'],
                                              catalog_magnitude=(o['Vmag'], o['e_Vmag']),
                                              cat_filter="V",
                                              id=o["recno"]))
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

    # TODO : change IDs
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
            object_list.append(SkyObject(SkyObject.star(ra=o['RAJ2000'],
                                                        dec=o['DEJ2000'],
                                                        catalog_magnitude=(o['Vmag'], o['e_Vmag']),
                                                        cat_filter="V",
                                                        id=i)))
            i += 1

        return object_list
