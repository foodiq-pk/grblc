from abc import ABC
from data_tier_placeholder.image import Image
from data_tier_placeholder.object import Object


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
        self.path = path

    def get_list(self):
        """returns list of objects based on selected handler and query"""
        raise NotImplementedError


class FileHandler(BasicHandler):
    """
    Handler for file and folder type queries.
    Only allows generation of lists of images
    """
    def __init__(self, path: str, query: str=None):
        super().__init__(query, path)

    def get_list(self):
        pass
    # TODO


class DatabaseHandler(BasicHandler):
    """
    Handler for interacting with database (sqlite3)
    Allows loading and saving images and object data from database file.
    """
    def __init__(self, path: str, query: str=None):
        super().__init__(query, path)

    def get_list(self):
        pass
    # TODO

    @staticmethod
    def save_objects_and_images(image_list: [Image], object_list: [Object]):
        """Saving to database file using SQL Alchemy"""
        pass
    # TODO

