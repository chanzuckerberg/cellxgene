"""
Drops and recreates all tables for local testing according to cellxgene_orm.py
"""
import glob
import os

from flask import current_app
from sqlalchemy import create_engine

from server.app.app import list_all_datasets
from server.db.cellxgene_orm import Base, CellxGeneDataset
from server.db.db_utils import DbUtils


def create_db(database_uri: str = "postgresql://postgres:test_pw@localhost:5432"):
    engine = create_engine(database_uri)
    print("Dropping tables")
    Base.metadata.drop_all(engine)
    print("Recreating tables")
    Base.metadata.create_all(engine)
    engine.connect().execute("""
            CREATE OR REPLACE RULE cxgdataset_table_ignore_duplicate_inserts AS
            ON INSERT TO cxgdataset
                WHERE EXISTS (
                  SELECT 1
                FROM cxgdataset
                WHERE name = NEW.name
            )
            DO INSTEAD NOTHING;
    """)
    print("Adding datasets")
    try:
        add_datasets_to_db(database_uri)
    except Exception as e:
        print(e)


def add_datasets_to_db(database_uri: str = "postgresql://postgres:test_pw@localhost:5432"):
    datasets_to_create = []
    # server_config = current_app.app_config.server_config
    # stored_datasets = list_all_datasets(server_config)
    stored_datasets = list_all_local_datasets("test/test_datasets/cxg-datasets")
    for dataset in stored_datasets:
        datasets_to_create.append(CellxGeneDataset(name=dataset))
    db = DbUtils(database_uri)

    db.session.add_all(datasets_to_create)
    db.session.commit()


def list_all_local_datasets(path_to_datasets):
    full_path_datasets = []
    files = glob.glob(f"{path_to_datasets}/*")
    for file in files:
        full_path_datasets.append(os.path.join(path_to_datasets, file))
    return full_path_datasets


if __name__ == "__main__":
    create_db()
