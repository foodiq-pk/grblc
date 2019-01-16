from sqlalchemy import Column, Integer, String, Float
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from config import Config

engine = create_engine(Config.DB_ENGINE)
Base = declarative_base(engine)


def init_session():
    Session = sessionmaker()
    session = Session(bind=engine)
    return session


class Frame(Base):
    __tablename__ = "images"

    id = Column(Integer, primary_key=True)
    time_jd = Column(Float)
    exposure = Column(Float)
    type = Column(String(10))
    path = Column(String)
    shift = Column(Float)
    shifterr = Column(Float)


class SObject(Base):
    __tablename__ = "skyobjects"

    id = Column(Integer, primary_key=True)
    ra = Column(Float)
    dec = Column(Float)
    catmag = Column(Float)
    catmagerr = Column(Float)


class Magnitude(Base):
    __tablename__ = "magnitudes"

    star_id = Column(Integer, primary_key=True)
    frame_id = Column(Integer, primary_key=True)
    mag = Column(Float, nullable=True)
    magerr = Column(Float, nullable=True)


# TODO optional?
class Shift(Base):
    __tablename__ = "shifts"

    star_id = Column(Integer, primary_key=True)
    frame_id = Column(Integer, primary_key=True)
    shift = Column(Float, nullable=True)
    shifterr = Column(Float, nullable=True)
