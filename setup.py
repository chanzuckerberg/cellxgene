from setuptools import setup, find_packages

with open("README.md", "rb") as fh:
    long_description = fh.read().decode()

with open("backend/server/requirements.txt") as fh:
    requirements = fh.read().splitlines()

with open("backend/server/requirements-prepare.txt") as fh:
    requirements_prepare = fh.read().splitlines()

setup(
    name="excellxgene",
    version="1.1.2",
    packages=find_packages(),
    url="https://github.com/czbiohub/cellxgene",
    license="MIT",
    author="Chan Zuckerberg Biohub",
    author_email="alexander.tarashansky@czbiohub.org",
    description="Web application for exploration of large scale scRNA-seq datasets, upgraded to enable end-to-end interactive analysis.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
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
    entry_points={"console_scripts": ["cellxgene = backend.server.cli.cli:cli"]},
    extras_require=dict(prepare=requirements_prepare),
)
