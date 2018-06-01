import json
from numpy import float32, integer
from flask import make_response, jsonify

class Float32JSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, float32):
            return float(obj)
        elif isinstance(obj, integer):
            return int(obj)
        return json.JSONEncoder.default(self, obj)

def make_payload(data, errormessage="", errorcode=200):
    """
    Creates JSON respons for requests
    :param data: json data
    :param errormessage: error message
    :param errorcode: http error code
    :return: flask json repsonse
    """
    error = False
    if errormessage:
        error = True
    # Questionable
    data = json.loads(json.dumps(data, cls=Float32JSONEncoder))
    return make_response(jsonify({
        "data": data,
        "status": {
            "error": error,
            "errormessage": errormessage,
        }
    }), errorcode)


