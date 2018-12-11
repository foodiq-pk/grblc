class Config:
    APERTURE = 4
    DARK_PATH = '/tmp/temp_master_dark.fits'
    FLAT_PATH = '/tmp/temp_master_flat.fits'
    FILE_DUMP = '/tmp/'   # Necessary for Pyraf photometry
    DB_ENGINE = "sqlite:///:memory:"
