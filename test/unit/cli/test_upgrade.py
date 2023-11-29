import unittest

from server.cli.upgrade import validate_version_str, split_version, version_gt


class CLIUpgradeTests(unittest.TestCase):
    """Test cases for CLI logic"""

    def test_validate_version_str(self):
        self.assertTrue(validate_version_str("0.1.2"))
        self.assertTrue(validate_version_str("0.1.2-RC", release_only=False))
        self.assertFalse(validate_version_str("0.1"))
        self.assertFalse(validate_version_str("0.1.2.3"))
        self.assertFalse(validate_version_str("0.1.2-RC"))

    def test_split_version_str(self):
        self.assertEqual(split_version("0.1.2"), [0, 1, 2])
        with self.assertRaises(AttributeError):
            split_version("0.1")

    def test_assert_verstion_gt(self):
        self.assertTrue(version_gt("1.0.0", "0.1.1"))
        self.assertTrue(version_gt("0.1.0", "0.0.1"))
        self.assertTrue(version_gt("0.0.1", "0.0.0"))
        self.assertFalse(version_gt("0.0.0", "0.0.0"))
        self.assertFalse(version_gt("0.0.0", "0.0.1"))
        self.assertFalse(version_gt("0.0.1", "0.1.0"))
        self.assertFalse(version_gt("0.1.1", "1.0.0"))
