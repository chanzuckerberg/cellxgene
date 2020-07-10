from sqlalchemy import (
    Column,
    create_engine,
    DateTime,
    ForeignKey,
    String,
    text,
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker


Base = declarative_base()

DEFAULT_DATETIME = text("now()")


class DBSessionMaker:
    def __init__(self, database_uri):
        self.engine = create_engine(database_uri, connect_args={"connect_timeout": 5})
        self.session_maker = sessionmaker(bind=self.engine)

    def session(self, **kwargs):
        return self.session_maker(**kwargs)


class CellxgGeneUser(Base):
    """
    A registered CellxGene user.
    Links a user to their annotations
    """

    __tablename__ = "user"

    id = Column(String, primary_key=True)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)

    # Relationships
    annotations = relationship("Annotation", back_populates="user")


class Annotation(Base):
    """
    An annotation is a link between a user, a dataset and tiledb dataframe. A user can have multiple annotations for a
    dataset, the most recent annotation (based on created_at) will be the default returned when queried
    """

    __tablename__ = "annotation"

    id = Column(String, primary_key=True)
    tiledb_uri = Column(String)
    user_id = Column(String, ForeignKey("user.id"), nullable=False)
    dataset_id = Column(String, ForeignKey("dataset.id"), nullable=False)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)

    # Relationships
    user = relationship("CellxgGeneUser", back_populates="annotations")
    dataset = relationship("Dataset", back_populates="annotations")


class Dataset(Base):
    """
    Datasets refer to cellxgene datasets stored in tiledb
    """

    __tablename__ = "dataset"

    id = Column(String, primary_key=True)
    name = Column(String)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    annotations = relationship("Annotation", back_populates="dataset")
