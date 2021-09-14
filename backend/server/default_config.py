import yaml

default_config = """
server:
  app:
    verbose: false
    debug: false
    host: localhost
    port : null
    open_browser: false
    force_https: false
    flask_secret_key: null
    generate_cache_control_headers: false

  single_dataset:
    # If datapath is set, then cellxgene with serve a single dataset located at datapath.
    datapath: null
    obs_names: null
    var_names: null
    about: null
    title: null

  data_locator:
    s3:
      # s3 region name.
      #   if true, then the s3 location is automatically determined from the datapath or dataroot.
      #   if false/null, then do not set.
      #   if a string, then use that value (e.g. us-east-1).
      region_name: true

  adaptor:
    anndata_adaptor:
      backed: false

  limits:
    column_request_max: 32
    diffexp_cellcount_max: null


dataset:
  app:
    # Scripts can be a list of either file names (string) or dicts containing keys src, integrity and crossorigin.
    # these will be injected into the index template as script tags with these attributes set.
    scripts: []
    # Inline scripts are a list of file names, where the contents of the file will be injected into the index.
    inline_scripts: []

  presentation:
    max_categories: 1000
    custom_colors: true

  user_annotations:
    enable: true
    type: local_file_csv
    local_file_csv:
      directory: null
      file: null            # annotations file name
      gene_sets_file: null  # gene sets file name
    gene_sets:
      readonly: false       # gene sets CRUD enabled/disabled

  embeddings:
    names : []

  diffexp:
    enable: true
    lfc_cutoff: 0.01
    top_n: 10

  X_approximate_distribution: auto

external:
  # You can retrieve configuration parameters from this config file, the environment,
  # the AWS secrets manager, or from the "cellxgene launch" command line arguments.
  # They are applied in that order, meaning that if a parameter is defined in more
  # than one location, the last one applied takes effect.

  # environment variables:
  # This section describes how to map environment variables to configuration parameters.
  # The format is a list defining an environment variable.
  # Each entry in the list is a dictionary with three entries:
  # name:  the name of the environment variable
  # path: the path within the cellxgene configuration to update.
  # required: (default=False) a boolean.  If true, then it is an error if the environment variable is not set.

  environment:
     - name: CXG_SECRET_KEY
       path: [server, app, flask_secret_key]
       required: false

  # AWS Secrets Manager
  # This section describes how to map aws secrets to configuration parameters.
  # The format is the region for the secrets manager, then a list of secrets.
  # each secret has a name, and a list of values.
  # Each entry in the list of values is a dictionary with three entries:
  # key:  the key of the aws secret.
  # path: the path within the cellxgene configuration to update.
  # required: (default=False) a boolean.  If true, then it is an error if the key does not exist in the secret.
  #
  # example:
  # aws_secrets_manager:
  #   region: us-west-2
  #    - name: my_first_secret
  #      values:
  #      - key: flask_secret_key
  #        path: [server, app, flask_secret_key]
  #        required: true
  #      - key: db_uri
  #        path: [dataset, user_annotations, db_uri]
  #        required: true

  aws_secrets_manager:
    region: null
    secrets: []
"""


def get_default_config():
    return yaml.load(default_config, Loader=yaml.Loader)
