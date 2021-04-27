import fastobo
import fsspec
import os

from flask import current_app, has_request_context

from backend.common.errors import OntologyLoadFailure, DisabledFeatureError
from backend.common.utils.type_conversion_utils import get_schema_type_hint_of_array
from backend.common.genesets import write_gene_sets_tidycsv, read_gene_sets_tidycsv, validate_gene_sets
from backend.common.utils.data_locator import DataLocator
from backend.common.utils.utils import path_join


class Annotations:
    """ baseclass for annotations, including ontologies and genesets """

    """ our default ontology is the PURL for the Cell Ontology.
    See http://www.obofoundry.org/ontology/cl.html """
    DefaultOnotology = "http://purl.obolibrary.org/obo/cl.obo"

    def __init__(self, config={}):
        self.ontology_data = None
        self.config = config

    def user_annotations_enabled(self):
        return self.config.get("user-annotations", False)

    def check_user_annotations_enabled(self):
        if not self.user_annotations_enabled():
            raise DisabledFeatureError("User annotations are disabled.")

    def load_ontology(self, path):
        """Load and parse ontologies - currently support OBO files only."""
        if path is None:
            path = self.DefaultOnotology

        try:
            with fsspec.open(path) as f:
                obo = fastobo.iter(f)
                terms = filter(lambda stanza: type(stanza) is fastobo.term.TermFrame, obo)
                names = [tag.name for term in terms for tag in term if type(tag) is fastobo.term.NameClause]
                self.ontology_data = names

        except FileNotFoundError as e:
            raise OntologyLoadFailure("Unable to find OBO ontology path") from e

        except SyntaxError as e:
            raise OntologyLoadFailure(f"{path}:{e.lineno}:{e.offset} OBO syntax error, unable to read ontology") from e

        except Exception as e:
            raise OntologyLoadFailure(f"{path}:Error loading OBO file") from e

    def get_schema(self, data_adaptor):
        schema = []
        labels = self.read_labels(data_adaptor)
        if labels is not None and not labels.empty:
            for col in labels.columns:
                col_schema = dict(name=col, writable=True)
                col_schema.update(get_schema_type_hint_of_array(labels[col]))
                schema.append(col_schema)

        return schema

    def set_collection(self, name):
        """set or create a new annotation collection"""
        raise NotImplementedError

    def read_labels(self, data_adaptor):
        """Return the labels as a pandas.DataFrame"""
        raise NotImplementedError

    def write_labels(self, df, data_adaptor):
        """Write the labels (df) to a persistent storage such that it can later be read"""
        raise NotImplementedError

    def update_parameters(self, parameters, data_adaptor):
        """Update configuration parameters that describe information about the annotations feature"""
        params = {}
        params["annotations_genesets_readonly"] = True
        params["annotations_genesets_name_is_read_only"] = True
        parameters.update(params)

    @staticmethod
    def gene_sets_to_csv(genesets):
        """
        Convert the internal genesets format (returned by read_gene_set) into
        the simple Tidy CSV.
        """
        from io import StringIO

        if isinstance(genesets, dict):
            genesets = genesets.values()

        with StringIO() as sio:
            write_gene_sets_tidycsv(sio, genesets)
            return sio.getvalue()

    @staticmethod
    def gene_sets_to_response(genesets):
        """
        Convert the internal genesets format (returned by read_gene_set) into
        the dict expected by the JSON REST API
        """
        return list(genesets.values())

    def read_gene_sets(self, data_adaptor, context=None):
        if has_request_context():
            if not current_app.auth.is_user_authenticated():
                return ({}, 0)

        gene_sets_uri_or_path = dataset_uri_to_geneset_uri(data_adaptor.data_locator.uri_or_path)

        server_config = data_adaptor.server_config
        region_name = None if server_config is None else server_config.data_locator__s3__region_name
        gene_sets_locator = DataLocator(gene_sets_uri_or_path, region_name=region_name)
        if not gene_sets_locator.exists():
            return ({}, 0)

        gene_sets = read_gene_sets_tidycsv(gene_sets_locator, context)
        schema = data_adaptor.get_schema()
        var_index = schema["annotations"]["var"].get("index", "index")
        var_names = set(data_adaptor.query_var_array(var_index))

        gene_sets = validate_gene_sets(gene_sets, var_names)
        return (gene_sets, 0)


def dataset_uri_to_geneset_uri(data_uri_or_path):
    """ given a dataset URI, return the associated gene set URI """
    data_basename = os.path.basename(data_uri_or_path)
    base, ext = os.path.splitext(data_basename)
    if ext is not None:  # strip extension, if any
        data_basename = base

    genesets_basename = f"{data_basename}-genesets.csv"
    gene_sets_uri_or_path = path_join(data_uri_or_path, "..", genesets_basename)

    return gene_sets_uri_or_path
