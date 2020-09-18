import unittest


class BaseConfigTest(unittest.TestCase):
    def test_mapping_creation(self):
        pass

    def test_attr_check__valid_type(self):
        pass
        # check set type and list of types

    def test_attr_check__invalid_type(self):
        pass
        # check set type and list of types

    def test_check_config_ensuers_all_values_type_checked(self):
        pass

    def test_config_update__raises_error_for_unknow_keys(self):
        pass

    def test_config_update__marks_updated_value_as_unvalidated(self):
        pass

    def test_changes_from_default_returns_list_of_none_default_config_values(self):
        pass
