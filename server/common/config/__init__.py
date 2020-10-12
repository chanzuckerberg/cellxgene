import logging
import os
import sys

from server.common.aws_secret_utils import get_secret_key
from server.common.data_locator import discover_s3_region_name

DEFAULT_SERVER_PORT = 5005
BIG_FILE_SIZE_THRESHOLD = 100 * 2 ** 20  # 100MB

