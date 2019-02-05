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

## Recipe

Follow these steps to create a release.

1.  Preparation:
    - python3.6 environment, and a cellxgene clone
    - Define the release version number, using [semantic versioning](https://semver.org/),
      and specifying all three digits (eg, 0.3.0)
    - Write the release title and release notes and add to
      [release notes document](https://docs.google.com/document/d/1KnHwkYfhyWO5H8BDcMu7y3ogjvq5Yi4OwpmZ8DB6w0Y/edit)
2.  Create a release branch, eg, `release-version`
3.  In the release branch:
    - Run `make release-stage-1 PART=[major | minor | patch]` where you choose major/minor/patch depending on which part of the version
      is being bumped (eg, 0.2.9->0.3 is minor).
4.  Commit and push the new branch
5.  Create a PR for the release.
    - [optional] As needed, conduct PR review.
6.  Merge to master
7.  Create Github release using the version number and release notes
    ([instructions](https://help.github.com/articles/creating-releases/)).
    - Draft new release
    - Type version name matching release version number from (1)
    - Select `master` as release branch (ensure you merged the release PR)
    - Type title `Release {version num}`
    - [optional] Check pre-release if this release is not ready for production
    - Publish Release
8. Publish to pypi by performing the following steps (assumes you that you have registered for pypi, 
and that you have write access to the cellxgene pypi package):
    - Build the distribution and upload to test pypi `make release-stage-2`
    - [optional] Test the test installation in a fresh virtual environment using `make install-release-test`
    - Upload the package to real pypi using `make release-stage-final`
    - [optional] Test the installation in a fresh virtual environment using
      `pip install cellxgene`
    - **Troubleshooting**:
        - Fails to upload to test.pypi: pypi doesn't allow you to reupload a release with the same version number, 
        if you accidentally burned a release number you want to use on prod, you have a couple options. 
            1) OPTION 1: Create distribution `make pydist`; test release locally `pip install dist/<release tarball>`;
             then upload to prod `make release-stage-final`. 
            2) OPTION 2: (DANGER) release directly to prod: `make release-burned`.  
            3) OPTION 3: If the release was burned on prod as well run from Step 3 again with option
            PART=patch until you get to an unburned version.
        - The release doesn't install or fails your tests when you install it: Delete it from pypi - Go to pypi.org, sign in, 
        go to the cellxgene package, click manage, then in the options drop down, click delete and 
        follow the instructions. You will not be able to use that release number again. If it is a minor bug 
        and not a major regression, you can just release a patch. 
    

The optional steps are for testing purposes, and are recommended
for publishing any major releases, and any releases that significantly
change the packaging (e.g. new bundled files, new dependencies, etc.)
