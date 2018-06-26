import json


def parse_schema(filename):
    with open(filename) as fh:
        schema = json.load(fh)
    return schema
