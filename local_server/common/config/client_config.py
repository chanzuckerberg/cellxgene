from local_server import display_version as cellxgene_display_version


def get_client_config(app_config, data_adaptor):
    """
    Return the configuration as required by the /config REST route
    """

    server_config = app_config.server_config
    dataset_config = data_adaptor.dataset_config
    annotation = dataset_config.user_annotations
    auth = server_config.auth

    # FIXME The current set of config is not consistently presented:
    # we have camalCase, hyphen-text, and underscore_text

    # make sure the configuration has been checked.
    app_config.check_config()

    # display_names
    title = app_config.get_title(data_adaptor)
    about = app_config.get_about(data_adaptor)

    display_names = dict(engine=data_adaptor.get_name(), dataset=title)

    # library_versions
    library_versions = {}
    library_versions.update(data_adaptor.get_library_versions())
    library_versions["cellxgene"] = cellxgene_display_version

    # links
    links = {"about-dataset": about}

    # parameters
    parameters = {
        "layout": dataset_config.embeddings__names,
        "max-category-items": dataset_config.presentation__max_categories,
        "obs_names": server_config.single_dataset__obs_names,
        "var_names": server_config.single_dataset__var_names,
        "diffexp_lfc_cutoff": dataset_config.diffexp__lfc_cutoff,
        "backed": server_config.adaptor__anndata_adaptor__backed,
        "disable-diffexp": not dataset_config.diffexp__enable,
        "enable-reembedding": dataset_config.embeddings__enable_reembedding,
        "annotations": False,
        "annotations_file": None,
        "annotations_dir": None,
        "annotations_genesets": True,  # feature flag
        "annotations_genesets_readonly": dataset_config.user_annotations__gene_sets__readonly,
        "annotations_genesets_summary_methods": ["mean"],
        "annotations_cell_ontology_enabled": False,
        "annotations_cell_ontology_obopath": None,
        "annotations_cell_ontology_terms": None,
        "custom_colors": dataset_config.presentation__custom_colors,
        "diffexp-may-be-slow": False,
    }

    # corpora dataset_props
    # TODO/Note: putting info from the dataset into the /config is not ideal.
    # However, it is definitely not part of /schema, and we do not have a top-level
    # route for data properties.  Consider creating one at some point.
    corpora_props = data_adaptor.get_corpora_props()
    if corpora_props and "default_embedding" in corpora_props:
        default_embedding = corpora_props["default_embedding"]
        if isinstance(default_embedding, str) and default_embedding.startswith("X_"):
            default_embedding = default_embedding[2:]  # drop X_ prefix
        if default_embedding in data_adaptor.get_embedding_names():
            parameters["default_embedding"] = default_embedding

    data_adaptor.update_parameters(parameters)
    if annotation:
        annotation.update_parameters(parameters, data_adaptor)

    # gather it all together
    client_config = {}
    config = client_config["config"] = {}
    config["displayNames"] = display_names
    config["library_versions"] = library_versions
    config["links"] = links
    config["parameters"] = parameters
    config["corpora_props"] = corpora_props
    config["limits"] = {
        "column_request_max": server_config.limits__column_request_max,
        "diffexp_cellcount_max": server_config.limits__diffexp_cellcount_max,
    }

    if dataset_config.app__authentication_enable and auth.is_valid_authentication_type():
        config["authentication"] = {
            "requires_client_login": auth.requires_client_login(),
        }
        if auth.requires_client_login():
            config["authentication"].update(
                {
                    # Todo why are these stored on the data_adaptor?
                    "login": auth.get_login_url(data_adaptor),
                    "logout": auth.get_logout_url(data_adaptor),
                }
            )

    return client_config


def get_client_userinfo(app_config, data_adaptor):
    """
    Return the userinfo as required by the /userinfo REST route
    """

    server_config = app_config.server_config
    dataset_config = data_adaptor.dataset_config
    auth = server_config.auth

    # make sure the configuration has been checked.
    app_config.check_config()

    if dataset_config.app__authentication_enable and auth.is_valid_authentication_type():
        userinfo = {}
        userinfo["userinfo"] = {
            "is_authenticated": auth.is_user_authenticated(),
            "username": auth.get_user_name(),
            "user_id": auth.get_user_id(),
            "email": auth.get_user_email(),
            "picture": auth.get_user_picture(),
        }
        return userinfo
