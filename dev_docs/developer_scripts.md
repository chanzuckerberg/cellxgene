# Developer convenience scripts

This document describes scripts for accelerating cellxgene development.

Paths are relative to the root project directory. If you need to know what
this is, run `PROJECT_ROOT=$(git rev-parse --show-toplevel); echo
$PROJECT_ROOT`.

## Project-level scripts

### Build

**Usage:** from the `$PROJECT_ROOT` directory run:

- `make build` builds whole app client and server
- `make build-client` runs webpack build
- `make build-for-server-dev` builds client and copies output directly into
  source tree (only for server devlopment)

### Clean

Deletes generated files.

**Usage:** from the `$PROJECT_ROOT` directory run:

- `make clean` cleans everything including node modules (means build with take
  a while
- `make clean-lite` cleans built directories
- `make clean-server` cleans source tree

### Distribution

Creates distribution for python module to upload to pypi.

**Usage:** from the `$PROJECT_ROOT` directory run:

- `make pydist` builds code and then builds sdist

### Release

See `release_process.md`.

### Development environment

Installs requirements files.

**Usage:** from the `$PROJECT_ROOT` directory run:

- `make dev-env` installs requirements and requirments-dev (for building code)

### Installing cellxgene packages

**Usage:** from the `$PROJECT_ROOT` directory:

- `install-dev` - installs from local source tree
- `install-release-test` - installs from test pypi
- `install-release` - installs from pypi
- `install-dist` - installs from local dist folder
- `uninstall` - uninstalls cellxgene

## Client-level scripts

### Running the client

#### start-frontend

**About** Serve the current client javascript independently from the `server` code.

**Requires**

- The server to be running. Best way to do this is with [backend_dev](#backend_dev).
- `make ci` to install the necessary node modules

**Usage:** from the `$PROJECT_ROOT/client` directory run `make start-frontend`

#### backend_dev

**About** This script enables FE developers to run the REST API necessary to
back the development server for the front end. It is intended to ensure that
the FE developer gets the current version of the backend with a single command
and no knowledge of python necessary. It creates and activates a virtual
environment and installs cellxgene from the current branch.

**Requires** `Python3.10+`, `virtual-env`, `pip`

**Usage:** from the `$PROJECT_ROOT` directory run `./scripts/backend_dev`

**Options:**

- In parallel, you can then launch the node development server to serve the
  current state of the FE with [`start-frontend`](#start-frontend), usually in
  a different terminal tab.
- You can also select a specific dataset using `DATASET=<dataset path> ./scripts/backend_dev`.
- You can also use `CXG_OPTIONS` to pass options to the `cellxgene launch`
  command, as in `CXG_OPTIONS='--disable-annotations' ./scripts/backend_dev`.

**Breakdown**

| command                                  | purpose                                                    |
| ---------------------------------------- | ---------------------------------------------------------- |
| python3.12 -m venv cellxgene             | creates cellxgene virtual environment                      |
| source cellxgene/bin/activate            | activates virtual environment                              |
| yes \| pip uninstall cellxgene \|\| true | uninstalls cellxgene (if installed)                        |
| pip install -e .                         | installs current local version of cellxgene                |
| cellxgene launch                         | launches cellxgene (must supply dataset as last parameter) |

### Client test scripts

Methods used to test the client javascript code

**Usage:** from the `$PROJECT_ROOT/client` directory run:

- `make unit-test` Runs all unit tests. It excludes any tests in the e2e
  folder. This is used by travis to run unit tests.
- `make smoke-test` Starts backend development server and runs end to end
  tests. This is what travis runs. It depends on the `e2e` and the
  `backend-dev` targets. One starts the server, the other runs the tests. If
  developing a front-end feature and just checking if tests pass, this is
  probabaly the one you want to run.
- `npm run e2e` Runs backend tests without starting the server. You will need to
  start the rest api separately with the pbmc3k.h5ad file. Note you can use
  the `JEST_ENV` environment variable to change how JEST runs in the browser.
  The test runs against `localhost:3000` by default. You can use the
  `CXG_URL_BASE` env variable to test non-localhost deployments of cellxgene.
