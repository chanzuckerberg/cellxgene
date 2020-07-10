import os
import string
import sys
from random import random

from numba.core.typing.builtins import Int
from sqlalchemy import func

from server.db.cellxgene_orm import CellxgGeneUser, Dataset, Annotation, Base
from server.db.create_db import create_db
from server.db.db_utils import DbUtils

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..", ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa


class TestDatabase:
    def __init__(self):
        create_db()
        self.db = DbUtils()
        self._populate_test_data()
        del self.db

    def _populate_test_data(self):
        self._create_test_user()
        self._create_test_dataset()
        self._create_test_annotation()

    def _create_test_user(self):
        user = CellxgGeneUser(id="test_user_id")
        self.db.session.add(user)
        self.db.session.commit()

    def _create_test_dataset(self):
        dataset = Dataset(
            id="test_dataset_id",
            name="test_dataset",
        )
        self.db.session.add(dataset)
        self.db.session.commit()

    def _create_test_annotation(self):
        annotation = Annotation(
            id="test_project_link_id",
            tiledb_uri="tiledb_uri",
            user="test_user_id",
            dataset="test_dataset_id"
        )
        self.db.session.add(annotation)
        self.db.session.commit()

    @staticmethod
    def get_random_string():
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(12))

    def _create_test_users(self, user_count: Int = 10):
        users = []
        for i in range(user_count):
            users.append(CellxgGeneUser(id=self.get_random_string))
        self.db.session.add_all(users)
        self.db.session.commit()


    def _create_test_datasets(self, dataset_count: Int = 10):
        datasets = []
        for i in range(dataset_count):
            datasets.append(Dataset(id=self.get_random_string, name=self.get_random_string()))
        self.db.session.add_all(datasets)
        self.db.session.commit()

    def _create_test_datasets(self, dataset_count: Int = 10):
        datasets = []
        for i in range(dataset_count):
            datasets.append(Dataset(id=self.get_random_string, name=self.get_random_string()))

        self.db.session.add_all(datasets)
        self.db.session.commit()

    def order_by_random(self, table: Base):
        return self.session.query(table).order_by(func.random()).first()

    def _create_test_annotations(self, annotation_count: Int = 10):
        annotations = []
        for i in range(annotation_count):
            dataset = self.order_by_random(Dataset)
            user = self.order_by_random(CellxgGeneUser)
            annotations.append(Annotation(
                id=self.get_random_string(),
                tiledb_uri=self.get_random_string(),
                user=user,
                dataset=dataset
            ))
        self.db.session.add_all(annotations)
        self.db.session.commit()

