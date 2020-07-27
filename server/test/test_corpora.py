
import unittest
import anndata
import json

from server.common.corpora import corpora_get_versions_from_anndata, corpora_is_version_supported, corpora_get_props_from_anndata
from server.test import PROJECT_ROOT


class CorporaAPITest(unittest.TestCase):
    def test_corpora_get_versions_from_anndata(self):
        adata = self._get_h5ad()

        if 'version' in adata.uns:
            del adata.uns['version']
        self.assertIsNone(corpora_get_versions_from_anndata(adata))

        # something bogus
        adata.uns['version'] = 99
        self.assertIsNone(corpora_get_versions_from_anndata(adata))

        # something legit
        adata.uns['version'] = {
            'corpora_schema_version': '0.0.0',
            'corpora_encoding_version': '9.9.9'
        }
        self.assertEqual(
            corpora_get_versions_from_anndata(adata),
            ['0.0.0', '9.9.9']
        )

    def test_corpora_is_version_supported(self):
        self.assertTrue(corpora_is_version_supported('1.0.0', '0.1.0'))
        self.assertFalse(corpora_is_version_supported('0.0.0', '0.1.0'))
        self.assertFalse(corpora_is_version_supported('1.0.0', '0.0.0'))

    def test_corpora_get_props_from_anndata(self):
        adata = self._get_h5ad()

        if 'version' in adata.uns:
            del adata.uns['version']
        self.assertIsNone(corpora_get_props_from_anndata(adata))

        # something bogus
        adata.uns['version'] = 99
        self.assertIsNone(corpora_get_props_from_anndata(adata))

        # unsupported version, but missing required values
        adata.uns['version'] = {
            'corpora_schema_version': '99.0.0',
            'corpora_encoding_version': '32.1.0'
        }
        with self.assertRaises(ValueError):
            corpora_get_props_from_anndata(adata)

        # legit version, but missing required values
        adata.uns['version'] = {
            'corpora_schema_version': '1.0.0',
            'corpora_encoding_version': '0.1.0'
        }
        with self.assertRaises(KeyError):
            corpora_get_props_from_anndata(adata)

        some_fields = {
            "version": {
                'corpora_schema_version': '1.0.0',
                'corpora_encoding_version': '0.1.0'
            },
            "title": "title",
            "layer_descriptions": "layer_descriptions",
            "organism": "organism",
            "organism_ontology_term_id": "organism_ontology_term_id",
            "project_name": "project_name",
            "project_description": "project_description",
            "contributors": json.dumps([{
                "contributors": "contributors"
            }]),
            "project_links": json.dumps([
                {
                    "link_name": "link_name",
                    "link_url": "link_url",
                    "link_type": "SUMMARY"
                }
            ])
        }
        for k in some_fields:
            adata.uns[k] = some_fields[k]
        some_fields["contributors"] = json.loads(some_fields["contributors"])
        some_fields["project_links"] = json.loads(some_fields["project_links"])
        self.assertEqual(corpora_get_props_from_anndata(adata), some_fields)

    def _get_h5ad(self):
        return anndata.read_h5ad(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
