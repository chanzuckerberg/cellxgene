import json
import shutil
import tempfile
import unittest
from http import HTTPStatus

import anndata
import requests

from backend.czi_hosted.common.corpora import (
    corpora_get_versions_from_anndata,
    corpora_is_version_supported,
    corpora_get_props_from_anndata,
)
from backend.test.test_czi_hosted.unit import start_test_server, stop_test_server
from backend.test import PROJECT_ROOT

VERSION = "v0.2"


class CorporaAPITest(unittest.TestCase):
    def test_corpora_get_versions_from_anndata(self):
        adata = self._get_h5ad()

        if "version" in adata.uns:
            del adata.uns["version"]
        self.assertIsNone(corpora_get_versions_from_anndata(adata))

        # something bogus
        adata.uns["version"] = 99
        self.assertIsNone(corpora_get_versions_from_anndata(adata))

        # something legit
        adata.uns["version"] = {"corpora_schema_version": "0.0.0", "corpora_encoding_version": "9.9.9"}
        self.assertEqual(corpora_get_versions_from_anndata(adata), ["0.0.0", "9.9.9"])

    def test_corpora_is_version_supported(self):
        self.assertTrue(corpora_is_version_supported("1.0.0", "0.1.0"))
        self.assertFalse(corpora_is_version_supported("0.0.0", "0.1.0"))
        self.assertFalse(corpora_is_version_supported("1.0.0", "0.0.0"))

    def test_corpora_get_props_from_anndata(self):
        adata = self._get_h5ad()

        if "version" in adata.uns:
            del adata.uns["version"]
        self.assertIsNone(corpora_get_props_from_anndata(adata))

        # something bogus
        adata.uns["version"] = 99
        self.assertIsNone(corpora_get_props_from_anndata(adata))

        # unsupported version, but missing required values
        adata.uns["version"] = {"corpora_schema_version": "99.0.0", "corpora_encoding_version": "32.1.0"}
        with self.assertRaises(ValueError):
            corpora_get_props_from_anndata(adata)

        # legit version, but missing required values
        adata.uns["version"] = {"corpora_schema_version": "1.0.0", "corpora_encoding_version": "0.1.0"}
        with self.assertRaises(KeyError):
            corpora_get_props_from_anndata(adata)

        some_fields = {
            "version": {"corpora_schema_version": "1.0.0", "corpora_encoding_version": "0.1.0"},
            "title": "title",
            "layer_descriptions": "layer_descriptions",
            "organism": "organism",
            "organism_ontology_term_id": "organism_ontology_term_id",
            "project_name": "project_name",
            "project_description": "project_description",
            "contributors": json.dumps([{"contributors": "contributors"}]),
            "project_links": json.dumps([{"link_name": "link_name", "link_url": "link_url", "link_type": "SUMMARY"}]),
        }
        for k in some_fields:
            adata.uns[k] = some_fields[k]
        some_fields["contributors"] = json.loads(some_fields["contributors"])
        some_fields["project_links"] = json.loads(some_fields["project_links"])
        self.assertEqual(corpora_get_props_from_anndata(adata), some_fields)

    def test_corpora_get_props_from_anndata_v110(self):
        adata = self._get_h5ad()

        if "version" in adata.uns:
            del adata.uns["version"]
        self.assertIsNone(corpora_get_props_from_anndata(adata))

        # legit version, but missing required values
        adata.uns["version"] = {"corpora_schema_version": "1.1.0", "corpora_encoding_version": "0.1.0"}
        with self.assertRaises(KeyError):
            corpora_get_props_from_anndata(adata)

        # Metadata following schema 1.1.0, which removes some fields relative to 1.1.0
        some_110_fields = {
            "version": {"corpora_schema_version": "1.0.0", "corpora_encoding_version": "0.1.0"},
            "title": "title",
            "layer_descriptions": "layer_descriptions",
            "organism": "organism",
            "organism_ontology_term_id": "organism_ontology_term_id",
        }
        for k in some_110_fields:
            adata.uns[k] = some_110_fields[k]
        self.assertEqual(corpora_get_props_from_anndata(adata), some_110_fields)

    def _get_h5ad(self):
        return anndata.read_h5ad(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")


class CorporaRESTAPITest(unittest.TestCase):
    """ Confirm endpoints reflect Corpora-specific features """

    @classmethod
    def setCorporaFields(cls, path):
        adata = anndata.read_h5ad(path)
        corpora_props = {
            "version": {"corpora_schema_version": "1.0.0", "corpora_encoding_version": "0.1.0"},
            "title": "PBMC3K",
            "contributors": json.dumps([{"name": "name"}]),
            "layer_descriptions": {"X": "raw counts"},
            "organism": "human",
            "organism_ontology_term_id": "unknown",
            "project_name": "test project",
            "project_description": "test description",
            "project_links": json.dumps(
                [{"link_name": "test link", "link_type": "SUMMARY", "link_url": "https://a.u.r.l/"}]
            ),
            "default_embedding": "X_tsne",
        }
        adata.uns.update(corpora_props)
        adata.write(path)

    @classmethod
    def setUpClass(cls):
        cls.tmp_dir = tempfile.TemporaryDirectory()
        src = f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad"
        dst = f"{cls.tmp_dir.name}/pbmc3k.h5ad"
        shutil.copyfile(src, dst)
        cls.setCorporaFields(dst)
        cls.ps, cls.server = start_test_server([dst])

    @classmethod
    def tearDownClass(cls):
        stop_test_server(cls.ps)
        cls.tmp_dir.cleanup()

    def setUp(self):
        self.session = requests.Session()
        self.url_base = f"{self.server}/api/{VERSION}/"

    def test_config(self):
        endpoint = "config"
        url = f"{self.url_base}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")

        result_data = result.json()
        self.assertIsInstance(result_data["config"]["corpora_props"], dict)
        self.assertIsInstance(result_data["config"]["parameters"], dict)

        corpora_props = result_data["config"]["corpora_props"]
        parameters = result_data["config"]["parameters"]

        self.assertEqual(corpora_props["version"]["corpora_schema_version"], "1.0.0")
        self.assertEqual(corpora_props["organism"], "human")
        self.assertEqual(parameters["default_embedding"], "tsne")
