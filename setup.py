from setuptools import setup, find_packages

with open("README.md", "rb") as fh:
    long_description = fh.read().decode()

with open("server/requirements.txt") as fh:
    requirements = fh.read().splitlines()

with open("server/requirements-prepare.txt") as fh:
    requirements_prepare = fh.read().splitlines()

with open("server/requirements-annotate.txt") as fh:
    requirements_annotate = fh.read().splitlines()

setup(
    name="cellxgene",
    version="1.1.0-rc.0",
    packages=find_packages(),
    url="https://github.com/chanzuckerberg/cellxgene",
    license="MIT",
    author="Chan Zuckerberg Initiative",
    author_email="cellxgene@chanzuckerberg.com",
    description="Web application for exploration of large scale scRNA-seq datasets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
    python_requires=">=3.6",
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
    extras_require=dict(prepare=requirements_prepare, annotate=requirements_annotate),
)
