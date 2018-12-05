from setuptools import setup, find_packages

with open("README.md", "rb") as fh:
    long_description = fh.read().decode()

with open("server/requirements.txt") as fh:
    requirements = fh.read().splitlines()

setup(
    name="cellxgene",
    version="0.3.0",
    packages=find_packages(),
    url="https://github.com/chanzuckerberg/cellxgene",
    license="MIT",
    author="Colin Megill, Charlotte Weaver",
    author_email="cweaver@chanzuckerberg.com",
    description="Web application for exploration of large scale scRNA-seq datasets",
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=requirements,
    include_package_data=True,
    zip_safe=False,
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ),
    entry_points={
        "console_scripts":
            ["cellxgene = server.cli.cli:cli"]
    },
    extras_require=dict(
        louvain=['python-igraph', 'louvain>=0.6'],
    ),
)
