# Code of conduct

We warmly welcome contributions from the community!

To ensure a welcoming experience for our entire community, this project adheres to the Contributor Covenant
[code of conduct](https://github.com/chanzuckerberg/.github/tree/master/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code. Please report unacceptable behavior
to opensource@chanzuckerberg.com.

If you have any questions about any of this stuff, just ask! :)

# Contributing ideas and issues

We'd love to hear from you! Please submit any bug reports and feature requests through [Github issues](https://github.com/chanzuckerberg/cellxgene/issues).

# Direct contributions

## Getting started

If you are interested in working on `cellxgene` development, you'll need to use git to make a copy of the [project repository](https://www.youtube.com/watch?v=A-4WltCTVms&list=PLe6EXFvnTV7-_41SpakZoTIYCgX4aMTdU&index=2&t=0s) and share your changes.
If you're new to git, we recommend [GitKraken](https://www.gitkraken.com/) for an intuitive interface.

We have several "rules" for contributions:

1. If your contribution is complex, adds new features, new UI design or otherwise warrants discussion, we highly recommend that you submit a github issue, and engage other contributors in a discussion about the details of your proposed PR. This will save you time in the long run, as many details and decisions can be hashed out ahead-of-time.

2. Please submit any direct contributions by [forking the repository](https://www.youtube.com/watch?v=Lb4yvfrX_7I&list=PLe6EXFvnTV7-_41SpakZoTIYCgX4aMTdU&index=3&t=9s), creating a feature branch, and [submitting a Pull Request](https://www.youtube.com/watch?v=2VX1ISk9XH8&list=PLe6EXFvnTV7-_41SpakZoTIYCgX4aMTdU&index=9&t=0s). **Do not submit pull requests against master.**

First, you'll need the following installed on your machine

- python 3.6+
- node and npm (we recommend using [nvm](https://github.com/creationix/nvm) if this is your first time with node)

Then clone the project

```
git clone https://github.com/chanzuckerberg/cellxgene.git
```

This is enough to get you started with editing documentation. If you'd like to contribute code:

Build the client web assets by calling `make` from inside the `cellxgene` folder

```
make
```

Install all requirements (we recommend doing this inside a [virtual environment](install))

```
pip install -e .
```

You can start the app while developing either by calling `cellxgene` or by calling `python -m server`. We recommend using the `--debug` flag to see more output, which you can include when reporting bugs.

If you have any questions about developing or contributing, come hang out with us by joining the [CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and posting in the `#cellxgene-dev` channel.

## Contributing code

This project has made a few key design choices:

- The front-end is built with [`regl`](https://github.com/regl-project/regl) (a webgl library), [`react`](https://reactjs.org/), [`redux`](https://redux.js.org/), [`d3`](https://github.com/d3/d3), and [`blueprint`](https://blueprintjs.com/docs/#core) to handle rendering large numbers of cells with lots of complex interactivity
- The app is designed with a client-server model that can support a range of existing analysis packages for Python-based backend computational tasks (currently built for [scanpy](https://github.com/theislab/scanpy))
- The client uses fast cross-filtering to handle selections and comparisons across subsets of data

Depending on your background and interests, you might want to contribute to the frontend, or backend, or both!

Please submit any direct contributions via a Pull Request. It'd be great for PRs to include test cases and documentation updates where relevant, though we know the core test suite is itself still a work in progress.

## Contributing documentation

The documentation is written in [markdown](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet), and lives in the directory `cellxgene/docs/posts`. You can directly edit or add to these files and submit a Pull Request as described above.

To preview your changes on your local machine, you'll need to install Jekyll and Ruby using [these instructions](https://jekyllrb.com/docs/installation/) (you don't have to know how to program in Ruby, just install it).

You can then preview your changes by running `cellxgene/docs$ bundle exec jekyll serve` and navigating to the url indicated in the terminal.
