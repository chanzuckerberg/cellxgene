import unittest

from parameterized import parameterized_class

from server.common.errors import DatasetAccessError
from test import FIXTURES_ROOT
from test.unit import app_config


@parameterized_class(
    ("data_locator", "backed", "X_approximate_distribution"),
    [
        (f"{FIXTURES_ROOT}/pbmc3k_16.h5ad", True, "auto"),  # 16 bit conversion tests
    ],
)
class AdaptorLoadErrorTest(unittest.TestCase):
    def test_float16_backed_raises_err(self):
        with self.assertRaises(DatasetAccessError):
            config = app_config(
                self.data_locator,
                backed=self.backed,
                extra_dataset_config=dict(X_approximate_distribution=self.X_approximate_distribution),
            )
