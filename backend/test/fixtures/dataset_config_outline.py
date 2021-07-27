f"""
dataset:
  app:
    scripts: {scripts} #list of strs (filenames) or dicts containing keys
    inline_scripts: {inline_scripts} #list of strs (filenames)

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
      gene_sets_file: {local_file_csv_gene_sets_file}
    gene_sets:
      readonly: {gene_sets_readonly}

  embeddings:
    names: {embedding_names}

  diffexp:
    enable: {enable_difexp}
    lfc_cutoff: {lfc_cutoff}
    top_n: {top_n}

  X_approximate_distribution: {X_approximate_distribution}
"""
