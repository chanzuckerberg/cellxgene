---
layout: default
title: Troubleshooting
description: Troubleshooting
---
# Troubleshooting tips & tricks

#### I tried to `pip install cellxgene` and got a weird error I don't understand

This may happen, especially as we work out bugs in our installation process! Please create a new [Github issue](https://github.com/chanzuckerberg/cellxgene/issues), explain what you did, and include all the error messages you saw. It'd also be super helpful if you call `pip freeze` and include the full output alongside your issue.

#### I have a BIG dataset - how can I make cellxgene run as fast as possible?

If your dataset requires gigabytes of disk space, you may need to select an appropriate storage format in order to effectively utilize `cellxgene`. Tips and tricks:

- `cellxgene` is optimized for columnar data access. For large datasets, format the expression matrix (`.X`) as either a [SciPy CSC sparse matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html) or a dense Numpy array (whichever creates a smaller `h5ad` file). If you are using `cellxgene prepare`, include the `--sparse` flag to ensure `.X` is formatted as a CSC sparse matrix (by default, `.X` will be a dense matrix).
- `cellxgene` start time is directly proportional to `h5ad` file size and the speed of your file system. Expect that large (eg, million cell) datasets will take minutes to load, even on relatively fast computers with a high performance local hard drive. Once loaded, exploring metadata should still be quick.
- If your dataset size exceeds the size of memory (RAM) on the host computer, differential expression calculations will be extremely slow (or fail, if you run out of virtual memory). In this case, we recommend running with the `--disable-diffexp` flag.

#### I'm following the developer instructions and get an error about "missing files and directories” when trying to build the client

This is likely because you do not have node and npm installed, we recommend using [nvm](https://github.com/creationix/nvm) if you're new to using these tools.
