import mlflow


class FakeModel(mlflow.pyfunc.PythonModel):
    def __init__(self, input_to_output: dict = {}):
        self.input_to_output = input_to_output

    def predict(self, model_input) -> None:
        # this stdout output is useful for validating the input in a test, noting that this model will be invoked in a
        # subprocess, so stdout is one means of communicating information back to the test code
        print(f"__MODEL_INPUT__={model_input.iloc[0][0]}")
