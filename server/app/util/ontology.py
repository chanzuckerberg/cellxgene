"""
Load and parse ontologies - currently support OBO files only.
"""
import fsspec
import fastobo
import traceback  # use built-in formatter for SyntaxError


""" our default ontology is the PURL for the Cell Ontology.  See http://www.obofoundry.org/ontology/cl.html """
DefaultOnotology = "http://purl.obolibrary.org/obo/cl.obo"


class OntologyLoadFailure(Exception):
    pass


def load_obo(path):
    """ given a URI or path, return an array of term names """
    if path is None:
        path = DefaultOnotology

    try:
        with fsspec.open(path) as f:
            obo = fastobo.iter(f)
            terms = filter(lambda stanza: type(stanza) is fastobo.term.TermFrame, obo)
            names = [tag.name for term in terms for tag in term if type(tag) is fastobo.term.NameClause]
            return names

    except FileNotFoundError as e:
        raise OntologyLoadFailure(f"Unable to find OBO ontology path: {path}") from e

    except SyntaxError as e:
        msg = ''.join(traceback.format_exception_only(SyntaxError, e))
        raise OntologyLoadFailure(msg) from e

    except Exception as e:
        raise OntologyLoadFailure(f"Error loading OBO file {path}") from e
