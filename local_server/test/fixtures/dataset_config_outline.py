f"""
dataset:
  app:
    scripts: {scripts} #list of strs (filenames) or dicts containing keys
    inline_scripts: {inline_scripts} #list of strs (filenames)

    about_legal_tos: {about_legal_tos}
    about_legal_privacy: {about_legal_privacy}

    authentication_enable: {authentication_enable}

  presentation:
    max_categories: {max_categories}
    custom_colors: {custom_colors}

  user_annotations:
    enable: {enable_users_annotations}
    type: {annotation_type}
    local_file_csv:
      directory: {local_file_csv_directory}
      file: {local_file_csv_file}
    ontology:
      enable: {ontology_enabled}
      obo_location: {obo_location}

  embeddings:
    names: {embedding_names}
    enable_reembedding: {enable_reembedding}

  diffexp:
    enable: {enable_difexp}
    lfc_cutoff: {lfc_cutoff}
    top_n: {top_n}
"""
