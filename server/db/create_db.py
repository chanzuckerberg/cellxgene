"""
Drops and recreates all tables according to cellxgene_orm.py
"""

import os
import sys

from sqlalchemy import create_engine

from server.db.cellxgene_orm import Base

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa


def create_db():
    engine = create_engine('postgresql://postgres:test_pw@localhost:5432')
    print("Dropping tables")
    Base.metadata.drop_all(engine)
    print("Recreating tables")
    Base.metadata.create_all(engine)


if __name__ == "__main__":
    create_db()
