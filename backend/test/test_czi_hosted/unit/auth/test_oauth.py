import unittest
import random
import time
import base64
import json
import requests

from flask import Flask, jsonify, make_response, request, redirect
from multiprocessing import Process

import jose
from backend.czi_hosted.common.config.app_config import AppConfig
from backend.test import FIXTURES_ROOT

# This tests the oauth authentication type.
# This test starts a cellxgene server and a mock oauth server.
# API requests to login and logout and get the userinfo are made
# to the cellxgene server, which then sends requests to the mock
# oauth server.

# number of seconds that the oauth token is valid
from backend.test.test_czi_hosted.unit import BaseTest

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


# function to launch the mock oauth server
def launch_mock_oauth(mock_port):
    mock_oauth_app.run(port=mock_port)


class AuthTest(BaseTest):
    @classmethod
    def setUpClass(cls):
        # The port that the mock oauth server will listen on
        cls.mock_port = random.randint(10000, 12000)
        cls.dataset_dataroot = FIXTURES_ROOT
        cls.mock_oauth_process = Process(target=launch_mock_oauth, args=(cls.mock_port,))
        cls.mock_oauth_process.start()

        # Verify that the mock oauth server is ready (accepting requests) before starting the tests.

        # The following lines are polling until the mock server is ready.
        # The issue is we are starting a mock oauth server, then we are starting a cellxgene server,
        # which will start making requests to the mock oauth server.
        # So there is a race condition because the mock oauth server needs to be ready before it gets requests.
        # We check to see if it is ready, and if not we wait 1 second, then try again.
        # If it gets to 5 seconds, which is shouldn't, we assume something has gone wrong and fail the test.
        server_okay = False
        for _ in range(5):
            try:
                response = requests.get(f"http://localhost:{cls.mock_port}/.well-known/jwks.json")
                if response.status_code == 200:
                    server_okay = True
                    break
            except:  # noqa: E722
                pass

            # wait one second and try again
            time.sleep(1)

        assert server_okay

    @classmethod
    def tearDownClass(cls):
        cls.mock_oauth_process.terminate()

    def auth_flow(self, app_config, cookie_key=None):

        app_config.update_server_config(
            app__api_base_url="local",
            authentication__type="oauth",
            authentication__params_oauth__oauth_api_base_url=f"http://localhost:{self.mock_port}",
            authentication__params_oauth__client_id="mock_client_id",
            authentication__params_oauth__client_secret="mock_client_secret",
            authentication__params_oauth__jwt_decode_options={"verify_signature": False, "verify_iss": False},
        )

        app_config.update_server_config(multi_dataset__dataroot=self.dataset_dataroot)
        app_config.complete_config()

        server = self.create_app(app_config)
        server.testing = True
        session = server.test_client()

        # auth datasets
        config = json.loads(session.get("/d/pbmc3k.cxg/api/v0.2/config").data)
        userinfo = json.loads(session.get("/d/pbmc3k.cxg/api/v0.2/userinfo").data)

        self.assertFalse(userinfo["userinfo"]["is_authenticated"])
        self.assertIsNone(userinfo["userinfo"]["username"])
        self.assertTrue(config["config"]["authentication"]["requires_client_login"])
        self.assertTrue(config["config"]["parameters"]["annotations"])

        login_uri = config["config"]["authentication"]["login"]
        logout_uri = config["config"]["authentication"]["logout"]

        self.assertEqual(login_uri, "http://localhost:5005/login?dataset=d/pbmc3k.cxg/")
        self.assertEqual(logout_uri, "http://localhost:5005/logout?dataset=d/pbmc3k.cxg/")

        response = session.get(login_uri)
        # check that the login redirect worked

        self.assertEqual(response.status_code, 302)

        config = json.loads(session.get("/d/pbmc3k.cxg/api/v0.2/config").data)
        userinfo = json.loads(session.get("/d/pbmc3k.cxg/api/v0.2/userinfo").data)

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
            userinfo = json.loads(session.get(f"/d/pbmc3k.cxg/api/v0.2/userinfo").data)
            self.assertTrue(userinfo["userinfo"]["is_authenticated"])
            self.assertEqual(userinfo["userinfo"]["username"], "fake_user")

            cookie = session.cookies.get(cookie_key)
            token = json.loads(base64.b64decode(cookie))
            access_token_after = token.get("access_token")
            id_token_after = token.get("id_token")

            self.assertNotEqual(access_token_before, access_token_after)
            self.assertNotEqual(id_token_before, id_token_after)

            # invalid cookie is rejected
            session.cookies.set(cookie_key, "TEST_" + cookie)
            self.assertTrue(cookie_key in session.cookies)
            response = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo")
            # this is not an error, the invalid cookie is just ignored.
            self.assertEqual(response.status_code, 200)
            userinfo = json.loads(response.data)
            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])

            # invalid id_token is rejected
            test_token = token
            test_token["id_token"] = "TEST_" + id_token_after
            encoded_cookie = base64.b64encode(json.dumps(test_token).encode()).decode()
            session.cookies.set(cookie_key, encoded_cookie)
            response = session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo")
            # this is not an error, the invalid id_token is just ignored.
            self.assertEqual(response.status_code, 200)
            userinfo = json.loads(response.data)
            self.assertFalse(userinfo["userinfo"]["is_authenticated"])
            self.assertIsNone(userinfo["userinfo"]["username"])

        r = session.get(logout_uri)
        # check that the logout redirect worked

        self.assertEqual(r.history[0].status_code, 302)
        self.assertEqual(r.url, f"{server}/d/pbmc3k.cxg/")
        config = json.loads(session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/config").data)
        userinfo = json.loads(session.get(f"{server}/d/pbmc3k.cxg/api/v0.2/userinfo").data)
        self.assertFalse(userinfo["userinfo"]["is_authenticated"])
        self.assertIsNone(userinfo["userinfo"]["username"])
        self.assertTrue(config["config"]["parameters"]["annotations"])

    @unittest.skip("turn on when we utilizing auth in the explorer")
    def test_auth_oauth_session(self):
        # test with session cookies
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(
            authentication__params_oauth__session_cookie=True,
        )
        self.auth_flow(app_config)

    @unittest.skip("turn on when we utilizing auth in the explorer")
    def test_auth_oauth_cookie(self):
        # test with specified cookie
        app_config = AppConfig()
        app_config.update_server_config(app__flask_secret_key="secret")
        app_config.update_server_config(
            authentication__params_oauth__session_cookie=False,
            authentication__params_oauth__cookie=dict(key="test_cxguser", httponly=True, max_age=60),
        )

        self.auth_flow(app_config, "test_cxguser")
