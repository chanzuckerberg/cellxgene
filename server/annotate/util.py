from pathlib import Path

from server.common.utils.data_locator import DataLocator


def fetch_model(model_url, model_cache_dir=None) -> Path:
    model_url_locator = DataLocator(model_url, local_cache_dir=model_cache_dir)
    if not model_url_locator.exists():
        raise ValueError(f"model file '{model_url}' not found")

    return Path(model_url_locator.path)
