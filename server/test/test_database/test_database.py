import unittest
from server.db.cellxgene_orm import CellxGeneUser, CellxGeneDataset, Annotation
from server.db.db_utils import DbUtils
from server.test.fixtures.database import TestDatabase


class DatabaseTest(unittest.TestCase):
    db = DbUtils("postgresql://postgres:test_pw@localhost:5432")

    @classmethod
    def setUpClass(cls) -> None:
        TestDatabase()

    @classmethod
    def tearDownClass(cls) -> None:
        del cls.db

    def test_user_creation(self):
        one_user = self.db.get(table=CellxGeneUser, entity_id='test_user_id')
        self.assertEqual(one_user.id, 'test_user_id')
        user_count = self.db.session.query(CellxGeneUser).count()
        self.assertGreater(user_count, 10)

    def test_dataset_creation(self):
        one_dataset = self.db.query(table_args=[CellxGeneDataset], filter_args=[CellxGeneDataset.name == 'test_dataset'])
        self.assertEqual(one_dataset[0].name, 'test_dataset')
        dataset_count = self.db.session.query(CellxGeneDataset).count()
        self.assertGreater(dataset_count, 10)

    def test_annotation_creation(self):
        one_annotation = self.db.query(table_args=[Annotation], filter_args=[Annotation.tiledb_uri == 'tiledb_uri'])[0]
        self.assertEqual(one_annotation.tiledb_uri, 'tiledb_uri')
        annotation_count = self.db.session.query(Annotation).count()
        self.assertGreater(annotation_count, 10)
