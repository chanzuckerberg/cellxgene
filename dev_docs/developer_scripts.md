## Node scripts in package.json

####backend-dev
 
**About** This scripts enables FE developers to run the rest API necessary to back the development server for the front end. It is intended to ensure that the FE developer gets the current version of the backend with a single command and no knowledge of python necessary. It creates and activates a virtual environment and installs cellxgene from the current branch. 

**Requires** Python3.6 - `virtual-env`, `pip`

**Usage** `npm run backend-dev <dataset>` then you can launch the node development server to serve the current state of the FE `npm run start` 

**Tips** Developers will probably want to run this in parallel with the node dev server. You can either do this by running each in a separate termanial window or by running `backend-dev` in the background (add an ` &` at the end of the command to run in the background: `npm run backend-dev <dataset> &`).  

**Breakdown** 

| command | purpose |
| -------- |----------|
| python3.6 -m venv cellxgene  | creates cellxgene virtual environment |
| source cellxgene/bin/activate | activates virtual environment |
| yes \| pip uninstall cellxgene \|\| true | uninstalls cellxgene (if installed) |
| pip install -e .. | installs current local version of cellxgene |
| cellxgene launch | launches cellxgene (must supply dataset as last parameter |

### Test scripts

#### test

**About** Run test locally. In order for this command to succeed you will need to give it a specific unit test to run. It won't pass if you run all tests as may be expected. This is because the unit tests and end to end (e2e) tests require different testing environments.

**Usage** `npm run test <test file or pattern>`

#### unit-test

**About** Runs all unit tests. It excludes any tests in the e2e folder. This is used by travis to run unit tests.

**Usage** `npm run unit-test`

#### smoke-test

**About** Starts backend development server and runs end to end tests. This is what travis runs. It depends on the `e2e` and the `start-server-for-test` node scripts. One starts the server, the other runs the tests. If developing a front-end feature and just checking if tests pass, this is probabaly the one you want to run.    

**Requirements** Must set env variable for port. In terminal `export JEST_ENV='dev'` `export JEST_CXG_PORT=5000`

**Usage** `npm run smoke-test`

#### e2e

**About** Runs backend tests without starting the server. You will need to start the rest api separately with the pbmc3k.h5ad file.

**Requirements** Must set env variable for port. In terminal `export JEST_ENV='dev'` `export JEST_CXG_PORT=5000`

**Usage** `npm run e2e`

## makefile

### Build commands
builds source code

```
build - builds whole app client and server
build-server - makes build dir and moves python and client files
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
