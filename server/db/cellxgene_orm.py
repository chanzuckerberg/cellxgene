import uuid
from datetime import datetime

from sqlalchemy import (
    Column,
    DateTime,
    ForeignKey,
    String,
)
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()


class CellxGeneUser(Base):
    """
    A registered CellxGene user.
    Links a user to their annotations
    """

    __tablename__ = "cxguser"

    id = Column(String, primary_key=True)
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)
    updated_at = Column(DateTime, nullable=False, default=datetime.utcnow, onupdate=datetime.utcnow)

    # Relationships
    annotations = relationship("Annotation", back_populates="cxguser")


class Annotation(Base):
    """
    An annotation is a link between a user, a dataset and tiledb dataframe. A user can have multiple annotations for a
    dataset, the most recent annotation (based on created_at) will be the default returned when queried
    """

    __tablename__ = "annotation"

    id = Column(String, primary_key=True)
    tiledb_uri = Column(String)
    user_id = Column(String, ForeignKey("cxguser.id"), nullable=False)
    dataset_id = Column(UUID, ForeignKey("cxgdataset.id"), nullable=False)
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)

    # Relationships
    cxguser = relationship("CellxGeneUser", back_populates="annotations")
    dataset = relationship("CellxGeneDataset", back_populates="annotations")


class CellxGeneDataset(Base):
    """
    Datasets refer to cellxgene datasets stored by cellxgene
    """

    __tablename__ = "cxgdataset"

    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4, unique=True, nullable=False)
    name = Column(String, unique=True, index=True)
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)
    annotations = relationship("Annotation", back_populates="dataset")
