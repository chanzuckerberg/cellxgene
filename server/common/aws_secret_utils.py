import logging

import boto3
from flask import json

from server.common.errors import SecretKeyRetrievalError


def get_secret_key(region_name, secret_name):
    session = boto3.session.Session()
    client = session.client(service_name="secretsmanager", region_name=region_name)

    try:
        get_secret_value_response = client.get_secret_value(SecretId=secret_name)
        if "SecretString" in get_secret_value_response:
            var = get_secret_value_response["SecretString"]
            secret = json.loads(var)
            return secret
    except Exception as e:
        logging.critical(f"Caught exception during get_secret_key, {e}", exc_info=True)
        raise SecretKeyRetrievalError

    return None
