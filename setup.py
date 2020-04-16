from setuptools import setup, find_packages


def file_lines(filepath):
    with open(filepath) as fh:
        lines = fh.read().splitlines()
    return lines


with open("README.md", "rb") as fh:
    long_description = fh.read().decode()

setup(
    name="cellxgene",
    version="0.15.0",
    packages=find_packages(),
    url="https://github.com/chanzuckerberg/cellxgene",
    license="MIT",
    author="Chan Zuckerberg Initiative",
    author_email="cellxgene@chanzuckerberg.com",
    description="Web application for exploration of large scale scRNA-seq datasets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=file_lines("server/requirements.txt"),
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Framework :: Flask",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: JavaScript",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={"console_scripts": ["cellxgene = server.cli.cli:cli"]},
    extras_require=dict(
        prepare=file_lines("server/requirements-prepare.txt"),
        hosted=file_lines("server/requirements-hosted.txt")
    ),
)
