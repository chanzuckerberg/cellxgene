import os
import shutil
import unittest
import random
from unittest import mock
import yaml

from local_server.test import FIXTURES_ROOT


def mockenv(**envvars):
    return mock.patch.dict(os.environ, envvars)


class ConfigTests(unittest.TestCase):
    tmp_fixtures_directory = os.path.join(FIXTURES_ROOT, "tmp_dir")

    @classmethod
    def tearDownClass(cls) -> None:
        shutil.rmtree(cls.tmp_fixtures_directory)

    @classmethod
    def setUpClass(cls) -> None:
        os.makedirs(cls.tmp_fixtures_directory)

    def custom_server_config(
        self,
        verbose="false",
        debug="false",
        host="localhost",
        port="null",
        open_browser="false",
        force_https="false",
        flask_secret_key="secret",
        auth_type="session",
        insecure_test_environment="false",
        index="false",
        allowed_matrix_types=[],
        max_cached_datasets=5,
        timelimit_s=5,
        dataset_datapath="null",
        obs_names="null",
        var_names="null",
        about="null",
        title="null",
        data_locater_region_name="us-east-1",
        anndata_backed="false",
        column_request_max=32,
        diffexp_cellcount_max="null",
        config_file_name="server_config.yaml",
    ):
        configfile = os.path.join(self.tmp_fixtures_directory, config_file_name)
        server_config_outline_path = os.path.join(FIXTURES_ROOT, "server_config_outline.py")
        with open(server_config_outline_path, "r") as config_skeleton:
            config = config_skeleton.read()
            server_config = eval(config)
        with open(configfile, "w") as server_config_file:
            server_config_file.write(server_config)
        return configfile

    def custom_app_config(
        self,
        verbose="false",
        debug="false",
        host="localhost",
        port="null",
        open_browser="false",
        force_https="false",
        flask_secret_key="secret",
        auth_type="session",
        index="false",
        allowed_matrix_types=[],
        max_cached_datasets=5,
        timelimit_s=5,
        dataset_datapath="null",
        obs_names="null",
        var_names="null",
        about="null",
        title="null",
        data_locater_region_name="us-east-1",
        anndata_backed="false",
        column_request_max=32,
        diffexp_cellcount_max="null",
        scripts=[],
        inline_scripts=[],
        authentication_enable="true",
        max_categories=1000,
        custom_colors="true",
        enable_users_annotations="true",
        annotation_type="local_file_csv",
        db_uri="null",
        hosted_file_directory="null",
        local_file_csv_directory="null",
        local_file_csv_file="null",
        ontology_enabled="false",
        obo_location="null",
        embedding_names=[],
        enable_reembedding="false",
        enable_difexp="true",
        lfc_cutoff=0.01,
        top_n=10,
        environment=None,
        aws_secrets_manager_region=None,
        aws_secrets_manager_secrets=[],
        config_file_name="app_config.yml",
    ):
        random_num = random.randrange(999999)
        configfile = os.path.join(self.tmp_fixtures_directory, config_file_name)
        server_config = self.custom_server_config(
            verbose=verbose,
            debug=debug,
            host=host,
            port=port,
            open_browser=open_browser,
            force_https=force_https,
            flask_secret_key=flask_secret_key,
            auth_type=auth_type,
            index=index,
            allowed_matrix_types=allowed_matrix_types,
            max_cached_datasets=max_cached_datasets,
            timelimit_s=timelimit_s,
            dataset_datapath=dataset_datapath,
            obs_names=obs_names,
            var_names=var_names,
            about=about,
            title=title,
            data_locater_region_name=data_locater_region_name,
            anndata_backed=anndata_backed,
            column_request_max=column_request_max,
            diffexp_cellcount_max=diffexp_cellcount_max,
            config_file_name=f"temp_server_config_{random_num}.yml",
        )
        dataset_config = self.custom_dataset_config(
            scripts=scripts,
            inline_scripts=inline_scripts,
            authentication_enable=authentication_enable,
            max_categories=max_categories,
            custom_colors=custom_colors,
            enable_users_annotations=enable_users_annotations,
            annotation_type=annotation_type,
            db_uri=db_uri,
            hosted_file_directory=hosted_file_directory,
            local_file_csv_directory=local_file_csv_directory,
            local_file_csv_file=local_file_csv_file,
            ontology_enabled=ontology_enabled,
            obo_location=obo_location,
            embedding_names=embedding_names,
            enable_reembedding=enable_reembedding,
            enable_difexp=enable_difexp,
            lfc_cutoff=lfc_cutoff,
            top_n=top_n,
            config_file_name=f"temp_dataset_config_{random_num}.yml",
        )
        external_config = self.custom_external_config(
            environment=environment,
            aws_secrets_manager_region=aws_secrets_manager_region,
            aws_secrets_manager_secrets=aws_secrets_manager_secrets,
            config_file_name=f"temp_external_config_{random_num}.yml",
        )

        with open(configfile, "w") as app_config_file:
            app_config_file.write(open(server_config).read())
            app_config_file.write(open(dataset_config).read())
            app_config_file.write(open(external_config).read())

        return configfile

    def custom_dataset_config(
        self,
        scripts=[],
        inline_scripts=[],
        authentication_enable="true",
        max_categories=1000,
        custom_colors="true",
        enable_users_annotations="true",
        annotation_type="local_file_csv",
        db_uri="null",
        hosted_file_directory="null",
        local_file_csv_directory="null",
        local_file_csv_file="null",
        ontology_enabled="false",
        obo_location="null",
        embedding_names=[],
        enable_reembedding="false",
        enable_difexp="true",
        lfc_cutoff=0.01,
        top_n=10,
        config_file_name="dataset_config.yml",
    ):
        configfile = os.path.join(self.tmp_fixtures_directory, config_file_name)
        dataset_config_outline_path = os.path.join(FIXTURES_ROOT, "dataset_config_outline.py")
        with open(dataset_config_outline_path, "r") as config_skeleton:
            config = config_skeleton.read()
            dataset_config = eval(config)
        with open(configfile, "w") as dataset_config_file:
            dataset_config_file.write(dataset_config)

        return configfile

    def custom_external_config(
        self,
        environment=None,
        aws_secrets_manager_region=None,
        aws_secrets_manager_secrets=[],
        config_file_name="external_config.yaml",
    ):
        # set to the default if environment is None
        if environment is None:
            environment = [
                dict(name="CXG_SECRET_KEY", path=["server", "app", "flask_secret_key"], required=False),
            ]
        external_config = {
            "external": {
                "environment": environment,
                "aws_secrets_manager": {"region": aws_secrets_manager_region, "secrets": aws_secrets_manager_secrets},
            }
        }

        configfile = os.path.join(self.tmp_fixtures_directory, config_file_name)
        with open(configfile, "w") as external_config_file:
            yaml.dump(external_config, external_config_file)
        return configfile
