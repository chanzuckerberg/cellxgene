from .mlflow_model_fixture import FakeModel


def _load_pyfunc(data_path):
    return FakeModel()
