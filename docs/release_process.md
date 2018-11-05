# cellxgene Release Process

_This document defines the release process for cellxgene_

## Overview

The goal of the release process is to publish an installable package
to PyPi, with a matching tagged release on github.

The release process should result in the following side-effects:

- Version number bump
- JS assets built & packaged, committed to the repo
- Tagged github release
- Publication to PyPi

## Process

Follow these steps to create a release.

1.  Preparation:
    - Define the release version number
    - Write the release title and release notes
2.  Create a release branch, eg, `release-version`
3.  In the release branch:
    - set the version number in `client/package.json`
    - set the version number in `setup.py`
    - build the JS asserts using `bin/build-client`
4.  Commit and push the new branch
5.  Create a PR for the release. _As needed, conduct PR review._
6.  Create Github release using the version number and release notes ([instructions](https://help.github.com/articles/creating-releases/)).
7.  Publish to PyPi
