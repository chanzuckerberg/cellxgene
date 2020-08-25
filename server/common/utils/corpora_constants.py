class CorporaConstants(object):
    REQUIRED_SIMPLE_METADATA_FIELDS = [
        "version",
        "title",
        "layer_descriptions",
        "organism",
        "organism_ontology_term_id",
        "project_name",
        "project_description",
    ]

    # The Corpora specification requires some values encoded as JSON due to the inability of AnnData to store complex
    # types.
    REQUIRED_JSON_ENCODED_METADATA_FIELD = ["contributors", "project_links"]

    OPTIONAL_SIMPLE_METADATA_FIELDS = ["preprint_doi", "publication_doi", "default_embedding", "default_field", "tags"]
