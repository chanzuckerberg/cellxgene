# Developer guidelines

### Requirements
- npm
- Python 3.6+
- Chrome

[See dev section of README](../README.md)

**All instructions are expected to be run from the top level cellxgene directory unless otherwise specified.**

## Running test suite
Client and server tests run on Travis CI for every push, PR, and commit to master on github. End to end tests run nightly on master only. 

### Unit tests
Steps to run the all unit tests:
1. Start in the project root directory
1. `make dev-env`
1. `make unit-test`

### End to end tests

End to end tests use two env variables:
* `JEST_ENV` - environment to run end to end tests. Default `dev`
    * `prod` - run headless with no slowdown, chromium will not open.
    * `dev` - opens chromimum, runs tests with minimal slowdown, close on exit.
    * `debug` - opens chromium, runs tests with 100ms slowdown, dev tools open, chrome stays open on exit.
* `JEST_CXG_PORT` - port that end to end tests are being run on. Default `3000` (client hosted port).

On CI the end to end tests are run with `JEST_ENV` set to `prod` using the `smoke-test` make target.

To run end to end tests as they will be run on CI
1. cellxgene should be built and installed as [specified in server dev](#install)
2. `export JEST_ENV='prod'`
3. `export JEST_CXG_PORT=5000`
4. Run `npm run --prefix client/ smoke-test`

Run end to end tests interactively during development
1. cellxgene should be installed as [specified in client dev](#install-1)
2. Follow [launch](#launch-1) instructions for client dev with dataset `example-dataset/pbmc3k`
3. Run `make smoke-test`
4. To debug a failing test `export JEST_ENV='debug'` and re-run.

## Server dev
### Install
* Build the client and put static files in place: `make build-for-server-dev`
* Install from local files: `make install-dev`

### Launch
* `cellxgene launch [options] <datafile>`

### Reloading
If you install cellxgene using `make install-dev` the server will be restarted every time you make changes on the server code. If changes affects the client, the browser must be reloaded.

### Linter

We use `yapf` to auto-format and lint python.

To auto-format code run `make fmt`. To run lint checks on the code run `make lint`.

### Test
If you would like to run the server tests individually, follow the steps below
1. Install development requirements `make dev-env`
1. Run `make unit-test` in the `server` directory.

### Tips
* Install in a virtualenv
* May need to rebuild/reinstall when you make client changes

## Client dev
### Install
1. Install prereqs for client: `make dev-env`
2. Install cellxgene server: `make install-dev` Caveat: this will not build the production client package - you must use the [server install](#install) instructions above to serve web assets.

### Launch
To launch with hot reloading you need to launch the server and the client separately. Node's hot reloading starts the client on its own node server and auto-refreshes when changes are made.
1. Launch server (the client relies on the REST API being available): `cellxgene launch [options] <datafile>`
2. Launch client: in `client/` directory run `npm run start`
3. Client will be served on localhost:3000

### Build
To build only the client: `make build-client`

### Linter
We use `eslint` to lint the code and `prettier` as our code formatter.

### Test

If you would like to run the client tests individually, follow the steps below in the `client` directory
1. For unit tests run `npm run unit-test` or `make unit-test`
1. For the smoke test run `npm run smoke-test` or `make smoke-test`

### Tips
* You can also install/launch the server side code from npm scrips (requires python3.6 with virtualenv) in `client/` directory run `npm run backend-dev`


