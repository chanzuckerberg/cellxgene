"""This is a simple script to ensure the custom requirements.txt do not violate
the server requirements.txt. A hosted cellxgene deployment may specify the exact
version requirements on all the modules, and may add additional modules.
This script is meant to aid in making that list of custom requirements easier to maintain.
If cellxgene adds a new dependency, or changes the version requirements of an existing
dependency, then this script can check if the custom requirements are still valid"""

import sys
import requirements
from packaging.version import Version
import pkg_resources


def check(expected, custom):
    """checks that the custom requirements meet all the requirements of the expected requirements.
    The custom set of requirements may contain additional entries than expected.
    The requirements in custom must all be exact (==).
    An expected requirement must be present in custom, and must match all the specs
    for that requirement.

    expected : name of the expected requirement.txt file
    custom : name of the custom requirements.txt file
    """
    edict = parse_requirements(expected)
    cdict = parse_requirements(custom)

    okay = True

    # cdict must only have exact requirements (==)
    for cname, cspecs in cdict.items():
        if len(cspecs) != 1 or cspecs[0][0] != "==":
            print(f"Error, spec must be an exact requirement {custom}: {cname} {str(cspecs)}")
            okay = False

    for ename, especs in edict.items():
        if ename not in cdict:
            print(f"Error, missing requirement from {custom}: {ename} {str(especs)}")
            okay = False
            continue

        cver = Version(cdict[ename][0][1])
        for espec in especs:
            rokay = check_version(cver, espec[0], Version(espec[1]))
            if not rokay:
                print(f"Error, failed requirement from {custom}: {ename} {espec}, {cver}")
                okay = False

    if okay:
        print("requirements check successful")
        sys.exit(0)
    else:
        sys.exit(1)


def parse_requirements(fname):
    """Read a requirements file and return a dict of modules name / specification"""
    try:
        with open(fname, "r") as fd:
            try:
                # pylint: disable=no-member
                rdict = {req.name: req.specs for req in requirements.parse(fd)}
            except pkg_resources.RequirementParseError:
                print(f"Unable to parse the requirements file: {fname}")
                sys.exit(1)
    except Exception as e:
        print(f"Unable to open file {fname}: {str(e)}")
        sys.exit(1)

    return rdict


# pylint: disable=too-many-return-statements
def check_version(cver, optype, ever):
    """
    Simple version check.
    Note: There is more complexity to comparing version (PEP440).
    However the use cases in cellxgene are limited, and do not require a general solution.
    """

    if optype == "==":
        return cver == ever
    if optype == "!=":
        return cver != ever
    if optype == ">=":
        return cver >= ever
    if optype == ">":
        return cver > ever
    if optype == "<=":
        return cver <= ever
    if optype == "<":
        return cver < ever

    print(f"Error, optype not handled: {optype}")
    return False


if __name__ == "__main__":
    check(sys.argv[1], sys.argv[2])
