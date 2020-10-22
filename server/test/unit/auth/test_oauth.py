import unittest
import random
import time
import base64
import json
import requests

from flask import Flask, jsonify, make_response, request, redirect
from multiprocessing import Process

import jose
from server.common.config.app_config import AppConfig
from server.test import FIXTURES_ROOT, test_server

# This tests the oauth authentication type.
# This test starts a cellxgene server and a mock oauth server.
# API requests to login and logout and get the userinfo are made
# to the cellxgene server, which then sends requests to the mock
# oauth server.

# number of seconds that the oauth token is valid
TOKEN_EXPIRES = 2

# Create a mocked out oauth token, which servers all the endpoints needed by the oauth type.
mock_oauth_app = Flask("mock_oauth_app")


@mock_oauth_app.route("/authorize")
def authorize():
    callback = request.args.get("redirect_uri")
    state = request.args.get("state")
    return redirect(callback + f"?code=fakecode&state={state}")


@mock_oauth_app.route("/oauth/token", methods=["POST"])
def token():
    now = time.time()
    expires_at = now + TOKEN_EXPIRES
    headers = dict(alg="RS256", kid="fake_kid")
    payload = dict(name="fake_user", sub="fake_id", email="fake_user@email.com", email_verified=True, exp=expires_at)
    jwt = jose.jwt.encode(claims=payload, key="mysecret", algorithm="HS256", headers=headers)
    r = {
        "access_token": f"access-{now}",
        "id_token": jwt,
        "refresh_token": f"random-{now}",
        "scope": "openid profile email",
        "expires_in": TOKEN_EXPIRES,
        "token_type": "Bearer",
        "expires_at": expires_at,
    }
    return make_response(jsonify(r))


@mock_oauth_app.route("/v2/logout")
def logout():
    return_to = request.args.get("returnTo")
    return redirect(return_to)


@mock_oauth_app.route("/.well-known/jwks.json")
def jwks():
    data = dict(
        alg="RS256",
        kty="RSA",
        use="sig",
        kid="fake_kid",
    )
    return make_response(jsonify(dict(keys=[data])))


# The port that the mock oauth server will listen on
PORT = random.randint(10000, 12000)


# function to launch the mock oauth server
def launch_mock_oauth():
    mock_oauth_app.run(port=PORT)


class AuthTest(unittest.TestCase):
    def setUp(self):
        self.dataset_dataroot = FIXTURES_ROOT
        self.mock_oauth_process = Process(target=launch_mock_oauth)
        self.mock_oauth_process.start()

    def tearDown(self):
        self.mock_oauth_process.terminate()

    def auth_flow(self, app_config, cookie_key=None):

        app_config.update_server_config(
            app__api_base_url="local",
            authentication__type="oauth",
            authentication__params_oauth__oauth_api_base_url=f"http://localhost:{PORT}",
            authentication__params_oauth__client_id="mock_client_id",
            authentication__params_oauth__client_secret="mock_client_secret",
            authentication__params_oauth__jwt_decode_options={"verify_signature": False, "verify_iss": False},
        )

        app_config.update_server_config(multi_dataset__dataroot=self.dataset_dataroot)
        app_config.complete_config()

        with test_server(app_config=app_config) as server:
            session = requests.Session()

            # auth datasets
            config = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config").json()
            userinfo = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo").json()

            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])
            self.assertTrue(config["config"]["authentication"]["requires_client_login"])
            self.assertTrue(config["config"]["parameters"]["annotations"])

            login_uri = config["config"]["authentication"]["login"]
            logout_uri = config["config"]["authentication"]["logout"]

            self.assertEqual(login_uri, f"{server}/login?dataset=d/pbmc3k.cxg/")
            self.assertEqual(logout_uri, f"{server}/logout?dataset=d/pbmc3k.cxg/")

            r = session.get(login_uri)
            # check that the login redirect worked
            self.assertEqual(r.history[0].status_code, 302)
            self.assertEqual(r.url, f"{server}/d/pbmc3k.cxg/")
            config = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config").json()
            userinfo = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo").json()
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["username"], "fake_user")
            self.assertEqual(userinfo["userinfo"]["email"], "fake_user@email.com")
            self.assertTrue(config["config"]["parameters"]["annotations"])

            if cookie_key:
                cookie = session.cookies.get(cookie_key)
                token = json.loads(base64.b64decode(cookie))
                access_token_before = token.get("access_token")
                id_token_before = token.get("id_token")

                # let the token expire
                time.sleep(TOKEN_EXPIRES + 1)

                # check that refresh works
                session.get(login_uri)
                userinfo = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo").json()
                self.assertTrue(userinfo["userinfo"]["is_authenticated"])
                self.assertEqual(userinfo["userinfo"]["username"], "fake_user")

                cookie = session.cookies.get(cookie_key)
                token = json.loads(base64.b64decode(cookie))
                access_token_after = token.get("access_token")
                id_token_after = token.get("id_token")

                self.assertNotEqual(access_token_before, access_token_after)
                self.assertNotEqual(id_token_before, id_token_after)

                # invalid cookie fails
                # (THIS CURRENTLY ERRORS ON THE SERVER, BUT SHOULD RETURN EMPTY USERINFO OR 401)
                session.cookies.set(cookie_key, "TEST_" + cookie)
                userinfo = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo").json()
                self.assertIsNone(userinfo.get("userinfo"))

                # invalid id_token fails
                test_token = token
                test_token["id_token"] = "TEST_" + id_token_before
                encoded_cookie = base64.b64encode(json.dumps(test_token))
                session.cookies.set(cookie_key, encoded_cookie)
                userinfo = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo").json()
                self.assertFalse(userinfo["userinfo"]["is_authenticated"])
                self.assertIsNone(userinfo["userinfo"]["username"])

            r = session.get(logout_uri)
            # check that the logout redirect worked
            self.assertEqual(r.history[0].status_code, 302)
            self.assertEqual(r.url, f"{server}/d/pbmc3k.cxg/")
            config = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config").json()
            userinfo = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo").json()
            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])
            self.assertTrue(config["config"]["parameters"]["annotations"])

    def test_auth_oauth_session(self):
        # test with session cookies
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(
            authentication__params_oauth__session_cookie=True,
        )
        self.auth_flow(app_config)

    def test_auth_oauth_cookie(self):
        # test with specified cookie
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(
            authentication__params_oauth__session_cookie=False,
            authentication__params_oauth__cookie=dict(key="test_cxguser", httponly=True, max_age=60),
        )

        self.auth_flow(app_config, "test_cxguser")
