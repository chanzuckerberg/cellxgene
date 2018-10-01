import json

from numpy import float32, integer


class Float32JSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, float32):
            return float(obj)
        elif isinstance(obj, integer):
            return int(obj)
        return json.JSONEncoder.default(self, obj)
