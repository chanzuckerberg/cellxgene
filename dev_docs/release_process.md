# cellxgene Release Process

_This document defines the release process for cellxgene_

## Overview

The goal of the release process is to publish an installable package
to PyPi, with a matching tagged release on github.

The release process should result in the following side-effects:

-   Version number bump, using [semantic versioning](https://semver.org/)
-   JS assets built & packaged, committed to the repo
-   Tagged github release
-   Publication to PyPi

Note all release tags pushed to GitHub MUST follow semantic versioning.

## Releasing a Major or Minor Version of cellxgene

Please scroll down the section below for how to release a patch version. Follow these steps to create a major or minor release.

1.  Preparation:
    -   python3.6 environment, and a cellxgene clone
    -   Define the release version number, using [semantic versioning](https://semver.org/), and specifying all three digits (e.g., 0.3.0)
    -   Write the release title and release notes and add to [release notes document](https://docs.google.com/document/d/1KnHwkYfhyWO5H8BDcMu7y3ogjvq5Yi4OwpmZ8DB6w0Y/edit)
2.  Create a release branch, eg, `release-version-0.16.0`
3.  In the release branch, run `make create-release-candidate PART=[major | minor | patch]` where you choose major/minor/patch depending on which part of the version is being bumped (e.g., `0.2.9` -> `0.3.0` is minor version bump). This will bump the version and create a release *candidate* version (i.e. `0.3.0-rc.0`).
4.  Commit and push the new branch. This will trigger tests to ensure that your branch isn't broken.
5. Upload the release candidate to Test PyPI by running the command `make release-candidate-to-test-pypi`. (Make sure you are registered for PyPI and Test PyPI and you have write access to the cellxgene PyPI package for both).
6. Verify the release candidate in a fresh virtual environment by running `make install-release-test` which installs the cellxgene build you just uploaded the Test PyPI.
7. If you find errors with the release candidate, run `make recreate-release-candidate` to increment the release candidate version (i.e. `0.3.0-rc.0` -> `0.3.0-rc.1`). Then go back to Steps 5 and 6 to re-upload and re-test the new release candidate.
8. If everything looks good, push the release to Test PyPI without the release candidate tag by running the command `make release-final-to-test-pypi` (i.e. `0.3.0-rc.1` -> `0.3.0`).
    - **NOTE:** Once you push the final release version to Test PyPI, you cannot ever re-upload the build again. If you need to make changes to the build, you will have to "burn" the version number and bump the part again and go back to step 1 with a brand new version number. For example, if you upload `0.3.0` to Test PyPI and realize there's a bug, you will have to create a new version `0.4.0` and there will be no `0.3.0` version of cellxgene. This is why testing the release candidate is very important.
9.  Create a PR for the release and conduct a PR review.
10.  Merge to the `main` branch.
11.  Publish to PyPI (prod) (assuming you that you have registered for PyPI, and that you have write access to the cellxgene pypi package) by running `make release-final`.
12. Test the installation in a fresh virtual environment by running `pip install --no-cache-dir cellxgene`.
13.  Create Github release using the version number and release notes ([instructions](https://help.github.com/articles/creating-releases/)):
     - Draft new release
     - Type version name matching release version number from (1)
     - Select `main` as release branch (ensure you merged the release PR)
     - Type title `Release {version num}`
     - [optional] Check pre-release if this release is not ready for production
     - Publish Release

The optional steps are for testing purposes, and are recommended for publishing any major releases, and any releases that significantly change the packaging (e.g. new bundled files, new dependencies, etc.)

### Releasing a Patch Version of cellxgene (special case)

To make a bugfix release (a point release/patch release) when there are already other changes in `main` we need to do a modified version of our release process. The difference is that instead of using `main` we are going make our release branch off of the tag for the release we want to patch. We cherrypick the commits that we want to include in the patch. Then instead of merging to `main`, we create the release directly off of the branch.

1.  (same as above) Preparation:
    -   python3.6 environment, and a cellxgene clone
    -   Define the release version number, using [semantic versioning](https://semver.org/), and specifying all three digits (e.g., 0.3.2) (for this you will update the last digit to represent a bugfix change).
    -   Write the release title and release notes and add to [release notes document](https://docs.google.com/document/d/1KnHwkYfhyWO5H8BDcMu7y3ogjvq5Yi4OwpmZ8DB6w0Y/edit)
2.  Create a release branch off of the tag for the release you want to update.
    -   Checkout the tag for the release you want to fix. For example, if we are fixing 0.9.0: `git checkout 0.9.0`.
    -   Create a branch from that tag. `git branch release-version-0.9.1`
3.  Cherrypick the commits that you want included in this patch.
    -   Test that the cherrypicked commits landed and fixed the issue locally.
    -   We **WILL NOT** merge this branch back into `main` as these commits should already exist in `main`.
4.  In the release branch (i.e. `release-version-0.9.1`), run `make create-release-candidate PART=patch` to bump the patch version and create the first release candidate (i.e. `0.9.1-rc.0`).
5. Run `make release-candidate-to-test-pypi` to upload the release candidate to Test PyPI.
6. Verify the release candidate in a fresh virtual environment by running `make install-release-test` which installs the cellxgene build you just uploaded the Test PyPI.
7. If you find errors with the release candidate, run `make recreate-release-candidate` to increment the release candidate version (i.e. `0.9.1-rc.0` -> `0.9.1-rc.1`). Then go back to Steps 5 and 6 to re-upload and re-test the new release candidate.
8. If everything looks good, push the final version of the release to Test PyPI without the release candidate tag by running the command `make release-final-to-test-pypi` (i.e. `0.9.1-rc.1` -> `0.9.1`).
    - **NOTE:** Once you push the final release version to Test PyPI, you cannot ever re-upload the build again. If you need to make changes to the build, you will have to "burn" the version number and bump the part again and go back to step 1 with a brand new version number. For example, if you upload `0.9.1` to Test PyPI and realize there's a bug, you will have to create a new version `0.9.2` and there will be no `0.9.1` version of cellxgene. This is why testing the release candidate is very important.
9.  Commit and push the new branch. DO NOT MAKE A PR OR MERGE TO `main`.
    -   Wait for release to pass the tests.
10.  Publish to PyPI (prod) (assuming you that you have registered for PyPI, and that you have write access to the cellxgene pypi package) by running `make release-final`.
11. Test the installation in a fresh virtual environment by running `pip install --no-cache-dir cellxgene`.
12.  Create Github release using the version number and release notes
    ([instructions](https://help.github.com/articles/creating-releases/)).
    -   Draft new release
    -   Type version name matching release version number from (1)
    -   [**_Different than above_**] Select the release-branch you pushed at step 5 as release branch
    -   Type title `Release {version num}`
    -   [optional] Check pre-release if this release is not ready for production
    -   Publish Release

## Troubleshooting

### Fails to upload to test.pypi

_PyPi doesn't allow you to reupload a release with the same version number_
If you accidentally burned a release number you want to use on prod, you have a few options:

1. OPTION 1: Create distribution `make pydist`; test release locally `pip install dist/<release tarball>`;
   then upload to prod `make release-final`.
2. OPTION 2: (DANGER) release directly to prod: `make release-directly-to-prod`.
3. OPTION 3: If the release was burned on prod as well run from Step 3 again with option PART=patch until you get to an unburned version.

### The release doesn't install or fails your tests when you install it

Delete it from pypi - Go to pypi.org -> sign in -> go to the cellxgene package -> click manage -> then in the options drop down click delete -> follow the instructions. You will not be able to use that release number again. If it is a minor bug and not a major regression, you can just release a patch.

### If you need to run the final upload to PyPI (prod) on a different computer than where you ran the command to upload to Test PyPI.

If you run `make release-final` without running `make release-final-to-test-pypi` first, the dist will not have been build on the computer running the final PyPI push. The solution is to run `make release-directly-to-prod`. This both builds the distribution files and then releases directly to prod pypi.org.

## Command Details

### Initial creation stage - `make create-release-candidate PART=[major | minor | patch]`

    1. Pip installs requirements-dev
    2. Bumps version by [PART] and creates the first release candidate.
    3. Deletes build directory, client/build, dist and cellxgene.egg-info
    4. Creates the package-lock.json

### Test PyPI upload stage - `make release-candidate-to-test-pypi`

    1. Pip installs requirements-dev
    2. Builds client and server
    3. Creates distribution release (sdist)
    4. Uploads to test.pypi.org
    
### Recreating release candidate stage(s) - `make recreate-release-candidate`

    1. Pip installs requirements-dev
    2. Bumps release candidate version number.
    3. Deletes build directory, client/build, dist and cellxgene.egg-info
    4. Creates the package-lock.json

### Penultimate stage, final release to Test PyPI - `make release-final-to-test-pypi`
    1. Pip installs requirements-dev
    2. Removes release candidate tag from the version number.
    3. Deletes build directory, client/build, dist and cellxgene.egg-info
    4. Creates the package-lock.json
    5. Pip installs requirements-dev
    6. Builds client and server
    7. Creates distribution release (sdist)
    8. Uploads to test.pypi.org

### Final stage - `make release-final`

    ** Does not build distribution **
    1. Uploads to pypi.org

### (DANGER) Release directly to prod `make release-directly-to-prod`

    ** builds distribution and uploads directly to prod **
    Only use this if you are directed to by the troubleshooting guide
    1. Pip installs requirements-dev
    2. Builds client and server
    3. Creates distribution release (sdist)
    4. Uploads to pypi.org
