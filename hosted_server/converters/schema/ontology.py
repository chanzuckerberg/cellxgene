"""Methods for working with ontologies and the OLS."""
from urllib.parse import quote_plus

import requests

OLS_API_ROOT = "http://www.ebi.ac.uk/ols/api"

# Curie means something like CL:0000001


def _ontology_name(curie):
    """Get the name of the ontology from the curie, CL or UBERON for example."""
    return curie.split(":")[0]


def _ontology_value(curie):
    """Get the id component of the curie, 0000001 from CL:0000001 for example."""
    return curie.split(":")[1]


def _double_encode(url):
    """Double url encode a url. This is required by the OLS API."""
    return quote_plus(quote_plus(url))


def _iri(curie):
    """Get the iri from a curie. This is a bit hopeful that they all map to purl.obolibrary.org"""
    if _ontology_name(curie) == "EFO":
        return f"http://www.ebi.ac.uk/efo/EFO_{_ontology_value(curie)}"
    return f"http://purl.obolibrary.org/obo/{_ontology_name(curie)}_{_ontology_value(curie)}"


class OntologyLookupError(Exception):
    """Exception for some problem with looking up ontology information."""


def _ontology_info_url(curie):
    """Get the to make a GET to to get information about an ontology term."""

    # If the curie is empty, just return an empty string. This happens when there is no
    # valid ontology value.
    if not curie:
        return ""
    else:
        return f"{OLS_API_ROOT}/ontologies/{_ontology_name(curie)}/terms/{_double_encode(_iri(curie))}"


def get_ontology_label(curie):
    """For a given curie like 'CL:1000413', get the label like 'endothelial cell of artery'"""

    url = _ontology_info_url(curie)

    if not url:
        return ""

    response = requests.get(url)

    if not response.ok:
        raise OntologyLookupError(
            f"Curie {curie} lookup failed, got status code {response.status_code}: {response.text}"
        )
    return response.json()["label"]


def lookup_candidate_term(label, ontology="cl", method="select"):
    """Lookup candidate terms for a label. This is useful when there is an existing label in a
    submitted dataset, and you want to find an appropriate ontology term.

    Args:
      label: the label to find ontology terms for
      ontology: the ontology to search in, cl or uberon or efo for example
      method: select or search. search provides much broader results

    Returns:
      list of (curie, label) tuples returned by OLS
    """
    # using OLS REST API [https://www.ebi.ac.uk/ols/docs/api]
    url = f"{OLS_API_ROOT}/{method}?q={quote_plus(label)}&ontology={ontology.lower()}"
    response = requests.get(url)

    if not response.ok:
        raise OntologyLookupError(
            f"Label {label} lookup failed, got status code {response.status_code}: {response.text}"
        )

    return [(r["obo_id"], r["label"]) for r in response.json()["response"]["docs"]]
