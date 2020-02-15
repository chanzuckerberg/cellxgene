# Developer convenience methods

## Makefile

Documentation for the `Makefile` targets in the project root directory.

### Build commands

builds source code

```
build - builds whole app client and server
build-cli - makes build dir and moves python and client files
build-client - runs webpack build
build-for-server-dev - builds client and copies output directly into source tree (only for server devlopment)
```

### Clean commands

deletes generated files

```
clean - cleans everything including node modules (means build with take a while
clean-lite - cleans built directories
clean-server - cleans source tree
```

### Dist commands

creates distribution for python module to upload to pypi

```
pydist - builds code and then builds sdist
```

### Release commands

see release_process.md

### Env commands

Installs requirements files

```
dev-env - installs requirements and requirments-dev (for building code)
gui-env - installs requirements and requirments-dev and requirments-gui (for building native app)
```

### install commands

Installs cellxgene from different locations

```
install - installs from local build directory
install-dev - installs from local source tree
install-release-test - installs from test pypi
install-release - installs from pypi
install-dist - installs from local dist folder
uninstall - uninstalls cellxgene
```

### gui

Commands for building the native app

```
build-assets - builds the image and icon assets for the gui to pull from
gui-spec-osx - creates the initial spec file for osx, do not run unless you are starting from scratch, one time only
gui-spec-windows - creates the initial spec file for windows, do not run unless you are starting from scratch, one time only
gui-build-osx - builds the app from the osx spec file
gui-build-windows - builds the app from the windows spec file
```
## Client Makefile

The following phony `make` targets in `client/Makefile` are convenience methods for getting you up and developing.

### Runner scripts

#### start-frontend

**About** Serve the current client javascript independently from the `server` code.

**Requires**
* The server to be running. Best way to do this is with `make backend-dev`
* `make ci` to install the necessary node modules

**Usage** `make start-frontend`

#### backend-dev

**About** This script enables FE developers to run the REST API necessary to back the development server for the front end. It is intended to ensure that the FE developer gets the current version of the backend with a single command and no knowledge of python necessary. It creates and activates a virtual environment and installs cellxgene from the current branch.

**Requires** Python3.6 - `virtual-env`, `pip`

**Usage** `make backend-dev`. Optionally, you can then launch the node development server to serve the current state of
the FE with `make start-frontend`. You can also select a specific dataset using `DATASET=<dataset path> make backend-dev`.
You can also use `CXG_OPTIONS` to pass options to the `cellxgene launch` command, as in
`CXG_OPTIONS='--experimental-annotations --experimental-annotations-file annotations.csv' make backend-dev`.

**Tips** Developers will probably want to run this in parallel with the node dev server. You can either do this by running each in a separate termanial window or by running `backend-dev` in the background (add an `&` at the end of the command to run in the background: `DATASET=<dataset> make backend-dev &`).

**Breakdown**

| command                                  | purpose                                                   |
| ---------------------------------------- | --------------------------------------------------------- |
| python3.6 -m venv cellxgene              | creates cellxgene virtual environment                     |
| source cellxgene/bin/activate            | activates virtual environment                             |
| yes \| pip uninstall cellxgene \|\| true | uninstalls cellxgene (if installed)                       |
| pip install -e .                        | installs current local version of cellxgene               |
| cellxgene launch                         | launches cellxgene (must supply dataset as last parameter) |

### Test scripts

#### test

**About** Run test locally. In order for this command to succeed you will need to give it a specific unit test to run. It won't pass if you run all tests as may be expected. This is because the unit tests and end to end (e2e) tests require different testing environments.

**Usage** `make test`

#### unit-test

**About** Runs all unit tests. It excludes any tests in the e2e folder. This is used by travis to run unit tests.

**Usage** `make unit-test`

#### smoke-test

**About** Starts backend development server and runs end to end tests. This is what travis runs. It depends on the `e2e` and the `backend-dev` targets. One starts the server, the other runs the tests. If developing a front-end feature and just checking if tests pass, this is probabaly the one you want to run.

**Usage** `make smoke-test`

#### e2e

**About** Runs backend tests without starting the server. You will need to start the rest api separately with the pbmc3k.h5ad file. Note you can use the `JEST_ENV` environment variable to change how JEST runs in the browser.

**Usage** `make e2e`
