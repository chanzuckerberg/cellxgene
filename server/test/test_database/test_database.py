import unittest

from server.db.cellxgene_orm import CellxgGeneUser, Dataset, Annotation
from server.db.create_db import create_db
from server.db.db_utils import DbUtils
from server.test.fixtures.database import TestDatabase


class AppConfigTest(unittest.TestCase):
    db = DbUtils()

    @classmethod
    def setUpClass(cls) -> None:
        TestDatabase()

    @classmethod
    def tearDownClass(cls) -> None:
        del cls.db

    def test_user_creation(self):
        one_user = self.db.get(table=CellxgGeneUser, entity_id='test_user_id')
        self.assertEqual(one_user.id, 'test_user_id')
        user_count = self.db.session.query(CellxgGeneUser).count()
        self.assertGreater(user_count, 10)

    def test_dataset_creation(self):
        one_dataset = self.db.get(table=Dataset, entity_id='test_dataset_id')
        self.assertEqual(one_dataset.id, 'test_dataset_id')
        dataset_count = self.db.session.query(Dataset).count()
        self.assertGreater(dataset_count, 10)

    def test_annotation_creation(self):
        one_annotation = self.db.get(table=Annotation, entity_id='test_annotation_id')
        self.assertEqual(one_annotation.id, 'test_annotation_id')
        annotation_count = self.db.session.query(Annotation).count()
        self.assertGreater(annotation_count, 10)
