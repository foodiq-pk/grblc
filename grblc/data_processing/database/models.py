from sqlalchemy import Column, String, Float

from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class Frame(Base):
    __tablename__ = "images"

    id = Column(String, primary_key=True)
    time_jd = Column(Float, nullable=False)
    exposure = Column(Float, nullable=False)
    type = Column(String(10), nullable=False)
    path = Column(String, nullable=False)
    additional = Column(String, nullable=True)
    # additional may contain info about src flux, sky, stacks, shift

    def __repr__(self):
        return f"Image('{self.id}','{self.time_jd:16.8f}','{self.exposure}','{self.type}')"


class SObject(Base):
    __tablename__ = "skyobjects"

    id = Column(String, primary_key=True)
    ra = Column(Float)
    dec = Column(Float)
    type = Column(String(10))
    additional = Column(String)
    # additional info may contain catalog magnitudes or julian date trigger or other related stuff

    def __repr__(self):
        return f"SkyObject('{self.id}','{self.ra}','{self.dec}','{self.type}','{self.additional}')"


