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

  authentication:
    # The authentication types may be "none", "session", "oauth"
    # none:  No authentication support, features like user_annotations must not be enabled.
    # session:  A session based userid is automatically generated. (no params needed)
    # oauth: oauth2 is used for authentication;  parameters are defined in params_oauth.
    type: session

    params_oauth:
       # url to the auth server
       api_base_url: null
       # client_id of this app
       client_id: null
       # the client_secret known to the auth server and this app
       client_secret: null
       # cellxgene server location;
       # the browser will be redirected to locations relative to this location during login and logout.
       # A value of None, indicates the client and server are on the localhost.  http://localhost:<port> will be used.
       callback_base_url: null

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
      file: null
    ontology:
      enable: false
      obo_location: null

  embeddings:
    names : []
    enable_reembedding: false

  diffexp:
    enable: true
    lfc_cutoff: 0.01
    top_n: 10

"""


def get_default_config():
    return yaml.load(default_config, Loader=yaml.Loader)
