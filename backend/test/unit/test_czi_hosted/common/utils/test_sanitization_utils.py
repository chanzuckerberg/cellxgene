import unittest

from backend.czi_hosted.common.utils.sanitization_utils import sanitize_values_in_list, sanitize_keys_in_dictionary


class TestSanitizationUtils(unittest.TestCase):
    def test__sanitize_values_in_list__not_strings_raises_exception(self):
        keys_to_sanitize = [1, 2, 3]

        with self.assertRaises(Exception) as exception_context:
            sanitize_values_in_list(keys_to_sanitize)

        self.assertIn("must contain all strings", str(exception_context.exception))

    def test__sanitize_values_in_list__not_all_strings_raises_exception(self):
        keys_to_sanitize = ["1", "2", 3]

        with self.assertRaises(Exception) as exception_context:
            sanitize_values_in_list(keys_to_sanitize)

        self.assertIn("must contain all strings", str(exception_context.exception))

    def test__sanitize_values_in_list__replace_non_ascii_character_with_underscore(self):
        keys_to_sanitize = ["abc.", "~abc", "a~b/c"]
        expected_sanitized_keys_dict = dict(zip(keys_to_sanitize, ["abc_", "_abc", "a_b_c"]))

        actual_sanitized_keys_dict = sanitize_values_in_list(keys_to_sanitize)

        self.assertEqual(expected_sanitized_keys_dict, actual_sanitized_keys_dict)

    def test__sanitize_keys_in_dictionary__replace_non_ascii_character_with_underscore(self):
        dictionary_to_sanitize = {"abc.": 3, "~abc": 4, "a~b/c": 5}
        expected_sanitized_dict = {"abc_": 3, "_abc": 4, "a_b_c": 5}

        actual_sanitized_dict = dictionary_to_sanitize
        sanitize_keys_in_dictionary(actual_sanitized_dict)

        self.assertEqual(expected_sanitized_dict, actual_sanitized_dict)

    def test__sanitize_keys_in_dictionary__non_string_key_raises_exception(self):
        dictionary_to_sanitize = {4: 3, "~abc": 4, "a~b/c": 5}

        with self.assertRaises(Exception) as exception_context:
            sanitize_keys_in_dictionary(dictionary_to_sanitize)

        self.assertIn("must contain all strings", str(exception_context.exception))

    def test__sanitize_keys_in_dictionary__replace_only_some_keys(self):
        dictionary_to_sanitize = {"abc": 3, "~abc": 4, "a~b/c": 5}
        expected_sanitized_dict = {"abc": 3, "_abc": 4, "a_b_c": 5}

        actual_sanitized_dict = dictionary_to_sanitize
        sanitize_keys_in_dictionary(actual_sanitized_dict)

        self.assertEqual(expected_sanitized_dict, actual_sanitized_dict)
