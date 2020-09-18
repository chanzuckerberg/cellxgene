import unittest


class TestDatasetConfig(unittest.TestCase):
    def test_init_datatset_config_sets_all_vars_from_default_config(self):
        # assert config error/key error not raised?
        pass

    def test_complete_config_handles_app_presentation_annotations_embeddings_and_diffexp(self):
        pass

    def test_handle_app_validates_attr_types_of_app_vars(self):
        pass

    def test_app_sets_script_vars(self):
        pass

    def test_app_raisies_error_when_src_script_not_passed_in(self):
        pass

    def test_handle_presentation_validates_attr_types_of_presentation_vars(self):
        pass

    def test_handle_user_annotations_validates_attr_type_of_annotation_vars(self):
        pass

    def test_handle_user_annotations_ensures_auth_is_enabled_with_valid_auth_type(self):
        pass

    def test_handle_user_annotations__adds_warning_message_if_annnotation_vars_set_when_annotations_disabled(self):
        pass

    def test_handle_user_annotations__instantiates_user_annotations_class_correctly(self):
        pass

    def test_handle_local_file_csv_annotations__checks_file_and_dir(self):
        pass

    def test_handle_local_file_csv_annotations__validates_passed_in_annotations_file(self):
        pass

    def test_handle_local_file_csv_annotations__validates_ontology_file(self):
        pass

    def test_handle_hosted_tiledb_annotations__validates_attr_type_of_annotation_vars(self):
        pass

    def test_handle_embeddings__validates_attr_types_of_embedding_vars(self):
        pass

    def test_handle_embeddings__instantiates_matrix_data_loader_and_checks_types(self):
        pass

    def test_handle_diffexp__validates_attr_types_of_diffexp_vars(self):
        pass

    def test_handle_diffexp__raises_warning_for_large_datasets(self):
        pass
