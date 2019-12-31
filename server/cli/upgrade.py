import click
import re

from github import Github, RateLimitExceededException
from requests.exceptions import ConnectionError
from .. import __version__

# Official SemVer regex: https://semver.org/
SEMVER_FORMAT = re.compile(
    r"^(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)"
    + r"(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)"
    + r"(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?"
    + r"(?:\+(?P<buildmetadata>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"
)


def log_upgrade_check():
    # Sanity-check that the CLI version is a properly-formatted SemVer string
    assert validate_version_str(__version__, release_only=False)

    # Get the current latest release
    try:
        github = Github(retry=0, timeout=5)
        cellxgene = github.get_organization("chanzuckerberg").get_repo("cellxgene")
        release_tag_generator = (release.tag_name for release in cellxgene.get_releases())
        latest_release = next(release_tag_generator, lambda ver_str: validate_version_str(ver_str))
        if version_gt(latest_release, __version__):
            click.echo(f"There's a new version of cellxgene available ({latest_release})!")
            click.echo("To upgrade, run the following: pip install --upgrade cellxgene\n")
    except (RateLimitExceededException, ConnectionError):
        click.echo("Upgrade check failed.\n")


def validate_version_str(version_str, release_only=True):
    """
    Test if a string conforms to SemVer format (https://semver.org/)
    :param version_str: a string to be validated
    :param release_only: only declare releases (not prereleases) valid
    :return: True if the version string is of a valid SemVer format else False
    """
    match = SEMVER_FORMAT.match(version_str)
    has_match = match is not None
    if release_only:
        return has_match and not match.group("prerelease")
    return has_match


def split_version(version_string):
    """
    Split a SemVer-formatted string into its component integers
    :param version_string: a SemVer string to be split
    :return: an array of three integers
    """
    match = SEMVER_FORMAT.match(version_string)
    return [int(match.group(group)) for group in ["major", "minor", "patch"]]


def version_gt(left_version, right_version):
    for left, right in zip(split_version(left_version), split_version(right_version)):
        if left > right:
            return True
        elif right > left:
            return False
    return False
