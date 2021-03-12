import string
import random


from sqlalchemy import func

from backend.server.db.cellxgene_orm import CellxGeneUser, CellxGeneDataset, Annotation, Base
from backend.server.db.create_db import create_db
from backend.server.db.db_utils import DbUtils


class TestDatabase:
    def __init__(self):
        local_db_uri = "postgresql://postgres:test_pw@localhost:5432"
        create_db(local_db_uri)
        self.db = DbUtils(local_db_uri)
        self._populate_test_data()
        self._populate_test_data_many()

    def _populate_test_data(self):
        self._create_test_user()
        self._create_test_dataset()
        self._create_test_annotation()

    def _populate_test_data_many(self):
        self._create_test_users()
        self._create_test_datasets()
        self._create_test_annotations()

    def _create_test_user(self):
        user = CellxGeneUser(id="test_user_id")
        user2 = CellxGeneUser(id="1234")
        self.db.session.add(user)
        self.db.session.add(user2)
        self.db.session.commit()

    def _create_test_dataset(self):
        dataset = CellxGeneDataset(name="test_dataset",)
        self.db.session.add(dataset)
        self.db.session.commit()

    def _create_test_annotation(self):
        dataset = self.db.query([CellxGeneDataset], [CellxGeneDataset.name == "test_dataset"],)[0]
        annotation = Annotation(tiledb_uri="tiledb_uri", user_id="test_user_id", dataset_id=str(dataset.id))
        self.db.session.add(annotation)
        self.db.session.commit()

    @staticmethod
    def get_random_string():
        letters = string.ascii_lowercase
        return "".join(random.choice(letters) for i in range(12))

    def _create_test_users(self, user_count: int = 10):
        users = []
        for i in range(user_count):
            users.append(CellxGeneUser(id=self.get_random_string()))
        self.db.session.add_all(users)
        self.db.session.commit()

    def _create_test_datasets(self, dataset_count: int = 10):
        datasets = []
        for i in range(dataset_count):
            datasets.append(CellxGeneDataset(name=self.get_random_string()))
        self.db.session.add_all(datasets)
        self.db.session.commit()

    def order_by_random(self, table: Base):
        return self.db.session.query(table).order_by(func.random()).first()

    def _create_test_annotations(self, annotation_count: int = 10):
        annotations = []
        for i in range(annotation_count):
            dataset = self.order_by_random(CellxGeneDataset)
            user = self.order_by_random(CellxGeneUser)
            annotations.append(
                Annotation(tiledb_uri=self.get_random_string(), user_id=user.id, dataset_id=str(dataset.id))
            )
        self.db.session.add_all(annotations)
        self.db.session.commit()
