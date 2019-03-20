class Config:
    APERTURE = 4
    DARK_PATH = '/tmp/temp_master_dark.fits'
    FLAT_PATH = '/tmp/temp_master_flat.fits'
    FILE_DUMP = '/tmp/'   # Necessary for Pyraf photometry
    DB_ENGINE = "sqlite:////tmp/test.db"

    # logging level

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