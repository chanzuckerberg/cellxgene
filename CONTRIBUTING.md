# contributing to cellxgene

We warmly welcome contributions from the community! Please submit any bug reports and feature requests through [Github issues](https://github.com/chanzuckerberg/cellxgene/issues). Please submit any direct contributions by forking the repository, creating a branch, and submitting a Pull Request. It'd be great for PRs to include test cases and documentation updates where relevant, though we know the core test suite is itself still a work in progress.

All code contributions and dependencies must be compatible with the project's [open-source license (MIT)](LICENSE.txt).

This project adheres to the Contributor Covenant
[code of conduct](https://github.com/chanzuckerberg/.github/tree/master/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code. Please report unacceptable behavior
to opensource@chanzuckerberg.com.

And finally, if you have any questions about any of this stuff, just ask! :)

## developer guide

This project has made a few key design choices

- The front-end is built with [`regl`](https://github.com/regl-project/regl) (a webgl library), [`react`](https://reactjs.org/), [`redux`](https://redux.js.org/), [`d3`](https://github.com/d3/d3), and [`blueprint`](https://blueprintjs.com/docs/#core) to handle rendering large numbers of cells with lots of complex interactivity
- The app is designed with a client-server model that can support a range of existing analysis packages for backend computational tasks (currently built for [scanpy](https://github.com/theislab/scanpy))
- The client uses fast cross-filtering to handle selections and comparisons across subsets of data

Depending on your background and interests, you might want to contribute to the frontend, or backend, or both!

If you are interested in working on `cellxgene` development, we recommend cloning the project from Gitub. First you'll need the following installed on your machine

- python 3.6+
- node and npm (we recommend using [nvm](https://github.com/creationix/nvm) if this is your first time with node)

Then clone the project

```
git clone https://github.com/chanzuckerberg/cellxgene.git
```

Build the client web assets by calling `make` from inside the `cellxgene` folder

```
make
```

Install all requirements (we recommend doing this inside a virtual environment)

```
pip install -e .
```

You can start the app while developing either by calling `cellxgene` or by calling `python -m server`. We recommend using the `--debug` flag to see more output, which you can include when reporting bugs.

If you have any questions about developing or contributing, come hang out with us by joining the [CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and posting in the `#cellxgene-dev` channel.
