import unittest
from backend.czi_hosted.db.cellxgene_orm import CellxGeneUser, CellxGeneDataset, Annotation
from backend.czi_hosted.db.db_utils import DbUtils
from backend.test.fixtures.database import TestDatabase


class DatabaseTest(unittest.TestCase):
    db = DbUtils("postgresql://postgres:test_pw@localhost:5432")

    @classmethod
    def setUpClass(cls) -> None:
        TestDatabase()

    @classmethod
    def tearDownClass(cls) -> None:
        del cls.db

    def test_user_creation(self):
        one_user = self.db.get(table=CellxGeneUser, entity_id="test_user_id")
        self.assertEqual(one_user.id, "test_user_id")
        user_count = self.db.session.query(CellxGeneUser).count()
        self.assertGreater(user_count, 10)

    def test_dataset_creation(self):
        one_dataset = self.db.query(
            table_args=[CellxGeneDataset], filter_args=[CellxGeneDataset.name == "test_dataset"]
        )
        self.assertEqual(one_dataset[0].name, "test_dataset")
        dataset_count = self.db.session.query(CellxGeneDataset).count()
        self.assertGreater(dataset_count, 10)

    def test_annotation_creation(self):
        one_annotation = self.db.query(table_args=[Annotation], filter_args=[Annotation.tiledb_uri == "tiledb_uri"])[0]
        self.assertEqual(one_annotation.tiledb_uri, "tiledb_uri")
        annotation_count = self.db.session.query(Annotation).count()
        self.assertGreater(annotation_count, 10)

    def test_get_most_recent_annotation_for_user_dataset(self):
        dataset_id = str(
            self.db.query(table_args=[CellxGeneDataset], filter_args=[CellxGeneDataset.name == "test_dataset"])[0].id
        )

        # have to commit separately because created_at time written on the db server
        self.db.session.add(Annotation(dataset_id=dataset_id, user_id="test_user_id", tiledb_uri="tiledb_uri_0"))
        self.db.session.commit()

        self.db.session.add(Annotation(dataset_id=dataset_id, user_id="test_user_id", tiledb_uri="tiledb_uri_1"))
        self.db.session.commit()

        self.db.session.add(Annotation(dataset_id=dataset_id, user_id="test_user_id", tiledb_uri="tiledb_uri_2"))
        self.db.session.commit()

        self.db.session.add(Annotation(dataset_id=dataset_id, user_id="test_user_id", tiledb_uri="tiledb_uri_3"))
        self.db.session.commit()

        self.db.session.add(Annotation(dataset_id=dataset_id, user_id="test_user_id", tiledb_uri="tiledb_uri_4"))
        self.db.session.commit()

        most_recent_annotation = self.db.query_for_most_recent(
            Annotation, [Annotation.dataset_id == dataset_id, Annotation.user_id == "test_user_id"]
        )

        self.assertEqual(most_recent_annotation.tiledb_uri, "tiledb_uri_4")
