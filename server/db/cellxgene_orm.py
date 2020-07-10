import enum
import os
import sys

from server.db.db_config import CellxGeneDbConfig
from sqlalchemy import (
    Boolean,
    Column,
    create_engine,
    DateTime,
    Enum,
    ForeignKey,
    ForeignKeyConstraint,
    Integer,
    String,
    text,
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

Base = declarative_base()
# deployment_stage = os.environ["DEPLOYMENT_STAGE"]

DEFAULT_DATETIME = text("now()")


class DBSessionMaker:
    def __init__(self):
        self.engine = create_engine(CellxGeneDbConfig().database_uri, connect_args={"connect_timeout": 5})
        self.session_maker = sessionmaker(bind=self.engine)

    def session(self, **kwargs):
        return self.session_maker(**kwargs)


class CellxgGeneUser(Base):
    """
    A registered Corpora user.
    Maintains user details such as contact information and access control settings.
    """

    __tablename__ = "user"

    id = Column(String, primary_key=True)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    updated_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)

    # Relationships
    annotations = relationship("annotation", back_populates="user")


class Annotation(Base):
    """
    A Corpora project represents an in progress or live submission of a lab experiment.
    DbProjects are associated with one or more single-cell datasets and links to external repositories.
    """

    __tablename__ = "annotation"

    id = Column(String, primary_key=True)
    tiledb_uri = Column(String)
    user = Column(ForeignKey("user.id"), nullable=False)
    dataset = Column(ForeignKey("dataset.id"), nullable=False)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)

    # Relationships
    user = relationship("CellxgGeneUser", back_populates="annotations")
    dataset = relationship("Dataset", back_populates="annotations")


class Dataset(Base):
    """
    Models a single experiment uploaded and processed by Corpora.
    Describes experiment metadata such as specimen and assay data.
    Related data files are represented by DbDataArtifacts.
    """

    __tablename__ = "dataset"

    id = Column(String, primary_key=True)
    name = Column(String)
    created_at = Column(DateTime, nullable=False, server_default=DEFAULT_DATETIME)
    annotations = relationship("annotation", back_populates="user")
