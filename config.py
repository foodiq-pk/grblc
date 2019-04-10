class Config:
    APERTURE = 4
    DARK_PATH = '/tmp'
    FLAT_PATH = '/tmp'
    FILE_DUMP = '/tmp/'   # Necessary for Pyraf photometry
    DB_ENGINE = "sqlite:////tmp/test.db"
    OVERWRITE = True
    DATA_DIR = "/home/foodiq/data/grbs/"

    # logging level
    @staticmethod
    def print():
        print("Aperture: ".format(Config.APERTURE))
        print("Dark path: ".format(Config.DARK_PATH))
        print("Flat path: ".format(Config.FLAT_PATH))
        print("Database: ".format(Config.DB_ENGINE))
        print("Overwrite: ".format(Config.OVERWRITE))
        print("Data directory: ".format(Config.DATA_DIR))

    @staticmethod
    def set_aperture(aperture: int):
        Config.APERTURE = aperture

    @staticmethod
    def set_dark_path(path: str):
        Config.DARK_PATH = path

    @staticmethod
    def set_flat_path(path: str):
        Config.FLAT_PATH = path

    @staticmethod
    def set_file_dump(path: str):
        Config.FILE_DUMP = path

    @staticmethod
    def set_db(path: str):
        Config.DB_ENGINE = path

    @staticmethod
    def set_data_dir(path: str):
        Config.DATA_DIR = path

    @staticmethod
    def set_overwrite(status: bool):
        Config.OVERWRITE = status
