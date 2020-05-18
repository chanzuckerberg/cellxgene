import yaml

default_config = """
# cellxgene configuration

server:
  verbose: false
  debug: false
  host: "127.0.0.1"
  port : null

  # Scripts can be a list of either file names (string) or dicts containing keys src, integrity and crossorigin.
  # these will be injected into the index template as script tags with these attributes set.
  scripts: []
  # Inline scripts are a list of file names, where the contents of the file will be injected into the index.
  inline_scripts: []

  open_browser: false
  about_legal_tos: null
  about_legal_privacy: null
  force_https: false
  flask_secret_key: null
  generate_cache_control_headers: false
  server_timing_headers: false
  csp_directives: null

presentation:
  max_categories: 1000
  custom_colors: true

multi_dataset:
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
  datapath: null
  obs_names: null
  var_names: null
  about: null
  title: null

user_annotations:
  enable: true
  type: local_file_csv
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

"""


def get_default_config():
    return yaml.load(default_config, Loader=yaml.Loader)
