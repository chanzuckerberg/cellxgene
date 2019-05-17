# Developer guidelines

### Requirements
- npm
- Python 3.6+
- Chrome

[See dev section of README](../README.md)


## Server Dev
### Install
* Build the client and put static files in place: `make build-for-server-dev`
* Install from local files: `make install-dev`

### Launch
* `cellxgene launch [options] <datafile>`

### Reloading
If you install cellxgene using `make install-dev` the server will be restarted every time you make changes on the server code. If changes affects the client, the browser must be reloaded.

### Linter
We use `flake8` to lint code. Travis CI runs `flake8 server`.

### Test
1. Install development requirements `pip install -r server/requirements-dev.txt`
2. Run tests `pytest server/test`

### Tips
* Install in a virtualenv
* May need to rebuild/reinstall when you make client changes

## Client Dev
### Install
1. Install prereqs for client: `npm install --prefix client/ client`
2. Install cellxgene server: `pip install -e .` Caveat: this will install necessary files for launching the cellxgene REST API but, you will need to build the client if you want to be able to view cellxgene in the browser when launched from the CLI.

### Launch
To launch with hot reloading you need to launch the server and the client separately. Node's hot reloading starts the client on it's own node server and auto-refreshes when changes are made.
1. Launch server (the client relies on the REST API being available): `cellxgene launch [options] <datafile>`
2. Launch client: `npm run start`
3. Client will be served on localhost:3000

### Build
To build only the client: `make build-client`

### Linter
We use `prettier` to lint the code.

### Test
Run `npm run unit-test`

### Tips
* You can also install/launch the server side code from npm scrips (requires python3.6 with virtualenv) `npm run backend-dev`

## Running tests
### Server
Install development requirements `pip install -r server/requirements-dev.txt`
Run tests `pytest server/test`

### Client
Run `npm run unit-test`

### Smoke tests

Smoke tests use two env variables
* `JEST_ENV` - environment to run smoke tests. Default `dev`
    * `prod` - run headless with no slowdown, chromium will not open.
    * `dev` - opens chromimum, runs tests with minimal slowdown, close on exit.
    * `debug` - opens chromium, runs tests with 100ms slowdown, dev tools open, chrome stays open on exit.
* `JEST_CXG_PORT` - port that smoke tests are being run on. Default `3000` (client hosted port).

Run smoke tests locally
1. cellxgene should be built and installed
2. run `npm run --prefix client/ smoke-test` from cellxgene directory

Run tests interactively
1. cellxgene should be installed
2. Follow [launch](#launch-2) instructions for client dev
3. Run `npm run --prefix client/ e2e`
4. To debug a failing test `export JEST_ENV='debug'` and re-run.



