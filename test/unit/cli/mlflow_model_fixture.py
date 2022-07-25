import shutil
from tempfile import TemporaryDirectory, mkstemp

import mlflow


def write_model(model) -> str:
    with TemporaryDirectory() as mlflow_model_dir:
        mlflow.pyfunc.save_model(mlflow_model_dir, python_model=model)
        return shutil.make_archive(mkstemp()[1], "zip", mlflow_model_dir)


class FakeModel(mlflow.pyfunc.PythonModel):
    def __init__(self, input_to_output: dict = {}):
        self.input_to_output = input_to_output

    def predict(self, context, model_input):
        # useful for validating the input in a test, noting that this model will be invoked in a subprocess
        print(f"MODEL_INPUT={model_input.iloc[0][0]}")

        return None  # [self.input_to_output[model_input[i]] for i in range(len(model_input))]
