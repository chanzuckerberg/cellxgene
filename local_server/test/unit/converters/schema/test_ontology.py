import json

import unittest
import unittest.mock

from server.converters.schema import ontology


class TestOntologyParsing(unittest.TestCase):
    def setUp(self):

        self.curies = ["UBERON:0002048", "HsapDv:0000174", "NCBITaxon:9606", "EFO:0008995"]

        self.names = ["UBERON", "HsapDv", "NCBITaxon", "EFO"]

        self.values = ["0002048", "0000174", "9606", "0008995"]

        self.iris = [
            "http://purl.obolibrary.org/obo/UBERON_0002048",
            "http://purl.obolibrary.org/obo/HsapDv_0000174",
            "http://purl.obolibrary.org/obo/NCBITaxon_9606",
            "http://www.ebi.ac.uk/efo/EFO_0008995",
        ]

        URL_ROOT = "http://www.ebi.ac.uk/ols/api/ontologies/"
        self.urls = [
            URL_ROOT + "UBERON/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FUBERON_0002048",
            URL_ROOT + "HsapDv/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FHsapDv_0000174",
            URL_ROOT + "NCBITaxon/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FNCBITaxon_9606",
            URL_ROOT + "EFO/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0008995",
        ]

        self.responses = {
            "UBERON:0002048": {
                "iri": "http://purl.obolibrary.org/obo/UBERON_0002048",
                "description": ["Respiration organ that develops as an outpocketing of the esophagus."],
                "label": "lung",
            },
            "HsapDv:0000174": {
                "iri": "http://purl.obolibrary.org/obo/HsapDv_0000174",
                "description": ["Infant stage that refers to an infant who is over 1 and under 2 months old."],
                "label": "1-month-old human stage",
            },
            "NCBITaxon:9606": {
                "iri": "http://purl.obolibrary.org/obo/NCBITaxon_9606",
                "description": None,
                "label": "Homo sapiens",
            },
            "EFO:0008995": {
                "iri": "http://www.ebi.ac.uk/efo/EFO_0008995",
                "description": [
                    (
                        '10X is a "synthetic long-read" technology and works by capturing a barcoded oligo-coated '
                        "gel-bead and 0.3x genome copies into a single emulsion droplet, processing the equivalent "
                        "of 1 million pipetting steps. Successive versions of the 10x chemistry use different "
                        "barcode locations to improve the sequencing yield and quality of 10x experiments."
                    )
                ],
                "label": "10X sequencing",
            },
        }

    def test_ontololgy_name(self):
        for curie, expected_name in zip(self.curies, self.names):
            self.assertEqual(ontology._ontology_name(curie), expected_name)

    def test_ontololgy_value(self):
        for curie, expected_value in zip(self.curies, self.values):
            self.assertEqual(ontology._ontology_value(curie), expected_value)

    def test_iri(self):
        for curie, expected_iri in zip(self.curies, self.iris):
            self.assertEqual(ontology._iri(curie), expected_iri)

    def test_ontology_info_url(self):
        for curie, expected_url in zip(self.curies, self.urls):
            self.assertEqual(ontology._ontology_info_url(curie), expected_url)

    def test_empty_ontology_info_url(self):
        self.assertEqual(ontology._ontology_info_url(""), "")


class TestOntologyLookup(unittest.TestCase):
    def setUp(self):
        self.responses = {
            "UBERON:0002048": {
                "iri": "http://purl.obolibrary.org/obo/UBERON_0002048",
                "description": ["Respiration organ that develops as an outpocketing of the esophagus."],
                "label": "lung",
            },
            "HsapDv:0000174": {
                "iri": "http://purl.obolibrary.org/obo/HsapDv_0000174",
                "description": ["Infant stage that refers to an infant who is over 1 and under 2 months old."],
                "label": "1-month-old human stage",
            },
            "NCBITaxon:9606": {
                "iri": "http://purl.obolibrary.org/obo/NCBITaxon_9606",
                "description": None,
                "label": "Homo sapiens",
            },
            "EFO:0008995": {
                "iri": "http://www.ebi.ac.uk/efo/EFO_0008995",
                "description": [
                    ('10X is a "synthetic long-read" technology and works by capturing a barcoded oligo-coated '
                     'gel-bead and 0.3x genome copies into a single emulsion droplet, processing the equivalent '
                     'of 1 million pipetting steps. Successive versions of the 10x chemistry use different barcode '
                     'locations to improve the sequencing yield and quality of 10x experiments.')
                ],
                "label": "10X sequencing",
            },
        }

        self.labels = {
            "UBERON:0002048": "lung",
            "HsapDv:0000174": "1-month-old human stage",
            "NCBITaxon:9606": "Homo sapiens",
            "EFO:0008995": "10X sequencing",
        }

    @unittest.mock.patch("requests.get")
    def test_lookup_label(self, mock_get):

        for curie, response in self.responses.items():
            mock_get.return_value.content = json.dumps(response)
            mock_get.return_value.json.return_value = response
            mock_get.return_value.status_code = 200

            label = ontology.get_ontology_label(curie)
            self.assertEqual(label, self.labels[curie])
