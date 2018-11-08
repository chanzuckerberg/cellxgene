# cellxgene release process

_This document defines the release process for cellxgene_

## Overview

The goal of the release process is to publish an installable package
to PyPi, with a matching tagged release on github.

The release process should result in the following side-effects:

- Version number bump, using semantic versioning
- JS assets built & packaged, committed to the repo
- Tagged github release
- Publication to PyPi

## Process

Follow these steps to create a release.

1.  Preparation:
    - Define the release version number, using [semantic versioning](https://semver.org/)
    - Write the release title and release notes
2.  Create a release branch, eg, `release-version`
3.  In the release branch:
    - set the version number in `client/package.json`
    - set the version number in `setup.py`
    - set the version number in `server/cli/cli.py`
    - build the JS asserts using `bin/build-client`
4.  Commit and push the new branch
5.  Create a PR for the release.
    - [optional] As needed, conduct PR review.
6.  Create Github release using the version number and release notes ([instructions](https://help.github.com/articles/creating-releases/)).
7.  Publish to pypi by performing the following steps
    (assumes you have `setuptools` and `twine` installed and that you have
    registered for pypi and have write access to the cellxgene pypi package)
    - build the distribution by calling
      `python setup.py sdist bdist_wheel`
      inside the top-level directory
    - [optional] upload the package to test pypi
      `twine upload --repository-url https://test.pypi.org/legacy/ dist/*`
    - [optional] test the test installation in a fresh virtual environment using
      `pip install --index-url https://test.pypi.org/simple/ cellxgene`
    - upload the package to real pypi using `twine upload dist/*`
    - [optional] test the installation in a fresh virtual environment using
      `pip install cellxgene`

The optional steps are for testing purposes, and are recommended
for publishing any major releases, and any releases that significantly
change the packaging (e.g. new bundled files, new dependencies, etc.)
