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
    server_timing_headers: false
    csp_directives: null

    # By default, cellxgene will serve api requests from the same base url as the webpage.
    # In general api_base_url and web_base_url will not need to be set.
    # There are two reasons to set these parameters:
    #  1. Oauth authentication is used; the oauth server will redirect back to the api_base_url after login,
    #     which then redirects back to the web_base_url.  If the web_base_url is not set, it will default to
    #     the api_base_url.  If oauth authentication is used, the api_base_url must be set.
    #     For a local test (where the server runs on "http://localhost:<port>"), then the api_base_url may be
    #     set to the string "local".
    #  2. The cellxgene deploymnent is in an environment where the webpage and api have
    #     different base urls.  In this case both api_base_url and web_base_url must be set.
    #     It is up to the server admin to ensure that the networking is setup correctly for this environment.
    api_base_url: null
    web_base_url: null

  authentication:
    # The authentication types may be "none", "session", "oauth"
    # none:  No authentication support, features like user_annotations must not be enabled.
    # session:  A session based userid is automatically generated. (no params needed)
    # oauth: oauth2 is used for authentication;  parameters are defined in params_oauth.
    type: session
    insecure_test_environment: false

    params_oauth:
       # url to the oauth server
       oauth_api_base_url: null
       # client_id of this app
       client_id: null
       # the client_secret known to the auth server and this app
       client_secret: null
       # jwt_decode_options, to specify non default decode options define
       # jwt_decode_options to be a dictionary with key/values described by
       # the options parameter of the jose.jwt.decode function:
       # (https://python-jose.readthedocs.io/en/latest/jwt/api.html)
       jwt_decode_options: null

       # if true, the jwt containing the id_token is stored in a session cookie
       session_cookie:  true

       # if session_cookie is false, then a regular cookie will be used.  In that case
       # the cookie will be defined by a dictionary of parameters.
       # The keys of the dictionary match the parameters of the flask set_cookie api
       # (https://flask.palletsprojects.com/en/1.1.x/api/), and with the same meaning.
       # legal keys:  key, max_age, expires, path, domain, secure, httponly, and samesite.
       cookie:  null

  multi_dataset:
    # If dataroot is set, then cellxgene may serve multiple datasets.  This parameter is not
    # compatible with single_dataset/datapath.
    # dataroot may be a string, representing the path to a directory or S3 prefix.  In this
    # case the datasets in that location are accessed from <server>/d/<datasetname>.
    # example:
    #     dataroot:  /path/to/datasets/
    # or
    #     dataroot:  s3://bucket/prefix/
    #
    # As an alternative, dataroot can be a dictionary, where a dataset key is associated with a base_url
    # and a dataroot.
    # example:
    #     dataroot:
    #         d1:
    #            base_url: set1
    #            dataroot: /path/to/set1_datasets/
    #         d2:
    #            base_url: set2/subdir
    #            dataroot: /path/to/set2_datasets/
    #
    # In this case, datasets can be accessed from <server>/set1/<datasetname> or
    # <server>/set2/subdir/<datasetname>.  It is possible to have different dataset configurations
    # for datasets accessed through different dataroots.  For example, in one dataroot, the
    # user annotations could be enabled, and in another dataroot they could be disabled.
    # To specify dataroot configurations, add a new top level dictionary to the config named
    # per_dataset_config. Within per_dataset_config create a dictionary for each dataroot to specialize
    # ("d1" or "d2" from the example).  Each of these dictionaries has the exact same form as the "dataset"
    # dictionary (see below).
    # When this approach is used, the values for each configuration option are checked in
    # this order: per_dataset_config/<key>, dataset, then the default values.
    #
    # example:
    #
    # per_dataset_config:
    #    d1:
    #       user_annotations:
    #           enable:  false
    #    d2:
    #       user_annotations:
    #           enable:  true

    dataroot: null

    # The index page when in multi-dataset mode:
    #   false or null:  this returns a 404 code
    #   true:  loads a test index page, which links to the datasets that are available in the dataroot
    #   string/URL:  redirect to this URL:  flask.redirect(config.multi_dataset__index)
    index: false

    # A list of allowed matrix types.  If an empty list, then all matrix types are allowed
    allowed_matrix_types: []

    matrix_cache:
      # The maximum number of datasets that may be opened at one time.  The least recently used dataset
      # is evicted from the cache first.
      max_datasets: 5

      # A matrix is automatically removed from the cache after timelimit_s number of seconds.
      # If timelimit_s is set to None, then there is no time limit.
      timelimit_s: 30

  single_dataset:
    # If datapath is set, then cellxgene with serve a single dataset located at datapath.  This parameter is not
    # compatible with multi_dataset/dataroot.
    datapath: null
    obs_names: null
    var_names: null
    about: null
    title: null

  diffexp:
    alg_cxg:
      # The number of threads to use is computed from: min(max_workers, cpu_multipler * cpu_count).
      # Where cpu_count is determined at runtime.
      max_workers: 64
      cpu_multiplier: 4

      # The target number of matrix elements that are evaluated
      # together in one thread.
      target_workunit: 16_000_000

  data_locator:
    s3:
      # s3 region name.
      #   if true, then the s3 location is automatically determined from the datapath or dataroot.
      #   if false/null, then do not set.
      #   if a string, then use that value (e.g. us-east-1).
      region_name: true

  adaptor:
    cxg_adaptor:
      # The key/values under tiledb_ctx will be used to initialize the tiledb Context.
      # If 'vfs.s3.region' is not set, then it will automatically use the setting from
      # data_locator / s3 / region_name.
      tiledb_ctx:
        sm.tile_cache_size:  8589934592
        sm.num_reader_threads:  32

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

    about_legal_tos: null
    about_legal_privacy: null

    # allow authentication support
    authentication_enable: true

  presentation:
    max_categories: 1000
    custom_colors: true

  user_annotations:
    enable: true
    type: local_file_csv
    hosted_tiledb_array:
        db_uri: null
        hosted_file_directory: null
    local_file_csv:
      directory: null
      file: null  # annotations file name
      genesets_file: null  # gene sets file name
    ontology:
      enable: false
      obo_location: null
    genesets:
      readonly: true

  embeddings:
    names : []
    enable_reembedding: false

  diffexp:
    enable: true
    lfc_cutoff: 0.01
    top_n: 10

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
     - name: CXG_OAUTH_CLIENT_SECRET
       path: [server, authentication, params_oauth, client_secret]
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
  #        path: [dataset, user_annotations, hosted_tiledb_array, db_uri]
  #        required: true
  #    - name:  my_auth_secret
  #      values:
  #      - key: client_secret
  #        path: [server, authentication, params_oauth, client_secret]
  #        required: true
  #      - key: client_id
  #        path: [server, authentication, params_oauth, client_id]
  #        required: true

  aws_secrets_manager:
    region: null
    secrets: []
"""


def get_default_config():
    return yaml.load(default_config, Loader=yaml.Loader)
