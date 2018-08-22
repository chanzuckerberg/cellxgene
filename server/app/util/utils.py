import json

from numpy import float32, integer
from flask import make_response, Response


class Float32JSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, float32):
            return float(obj)
        elif isinstance(obj, integer):
            return int(obj)
        return json.JSONEncoder.default(self, obj)


def make_payload(data, errorcode=200):
    """
    Creates JSON respons for requests
    :param data: json data
    :param errorcode: http error code
    :return: flask json repsonse
    """
    # Questionable
    data = json.loads(json.dumps(data, cls=Float32JSONEncoder))
    return make_response(data, errorcode)


def make_streaming_response(data_generator, errorcode=200, content_type="application/json"):
    # TODO headers
    return Response(data_generator, status=errorcode, content_type=content_type)
