f"""server:
  app:
    verbose: {verbose}
    debug: {debug}
    host: {host}
    port: {port}
    open_browser: {open_browser}
    force_https: {force_https}
    flask_secret_key: {flask_secret_key}
    generate_cache_control_headers: {generate_cache_control_headers}
    server_timing_headers: {server_timing_headers}
    csp_directives: {csp_directives}
    api_base_url: {api_base_url}
    web_base_url: {web_base_url}
  authentication:
    type: {auth_type}
    insecure_test_environment: {insecure_test_environment}
    params_oauth:
      oauth_api_base_url: {oauth_api_base_url}
      client_id: {client_id}
      client_secret: {client_secret}
      jwt_decode_options: {jwt_decode_options}
      session_cookie:  {session_cookie}
      cookie:  {cookie}

  multi_dataset:
    dataroot: {dataroot}
    index: {index}
    allowed_matrix_types: {allowed_matrix_types}
    matrix_cache:
      max_datasets: {max_cached_datasets}
      timelimit_s: {timelimit_s}

  single_dataset:
    datapath: {dataset_datapath}
    obs_names: {obs_names}
    var_names: {var_names}
    about: {about}
    title: {title}

  diffexp:
    alg_cxg: # number of threads to use is computed from: min(max_workers, cpu_multipler * cpu_count)
      max_workers: {diffexp_max_workers}
      cpu_multiplier: {cpu_multiplier}
      target_workunit: {target_workunit}  # The target number of matrix elements that are evaluated in one thread.

  data_locator:
    s3:
      region_name: {data_locater_region_name}

  adaptor:
    cxg_adaptor:
      tiledb_ctx:
        sm.tile_cache_size:  {cxg_tile_cache_size}
        sm.num_reader_threads:  {cxg_num_reader_threads}

    anndata_adaptor:
      backed: {anndata_backed}

  limits:
    column_request_max: {column_request_max}
    diffexp_cellcount_max: {diffexp_cellcount_max}
"""
