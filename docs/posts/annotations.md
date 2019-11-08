---
layout: default
title: annotations
description: Creating annotations
---

# Creating annotations in cellxgene

We are _piloting_ a new feature in cellxgene that enables users to create and edit categorical annotations within the app.

## Quick start for annotations  

You can enable this experimental feature like so:

`cellxgene launch mydata.h5ad --experimental-label-file myannotations.csv`

Any annotations you create in the application will be autosaved in `myannotations.csv`.

## Data management

To preserve data provenance, **`cellxgene` does not alter the input h5ad file**.

Rather, newly-created annotations are saved in the specified CSV file. If the specified file already exists, the previously-created annotations will be loaded as mutable (changeable) values and the CSV will be updated (overwritten) with any edits made. This allows you to create annotations and then come back later to continue working on them.

Once you're finished with your annotations, you should finalize and preserve your work by merging your `csv` into your main `h5ad` file.

You can do so like this:  
```
import pandas as pd
import scanpy as sc

new_annotations = pd.read_csv('myannotations.csv',
                         comment='#',
                         dtype='category',
                         index_col=0)
anndata = sc.read('mydata.h5ad')
anndata.obs = anndata.obs.join(new_annotations)
```

## FAQ

### How do I know my annotations are saved?

`cellxgene` autosaves any changes made to your annotations every 3 seconds.

### What about creating continuous annotations?  

Continuous metadata is important! However, these values (e.g., pseudotime) are the result of statistical analyses that are beyond cellxgene's visualization- and exploration-focused scope. We do, of course, provide visualization of continuous metadata values computed elsewhere and stored in `anndata.obs`.

### I keep getting weird index errors when trying to join my annotations to my anndata??  
This is most likely because the h5ad file you are working with is not the original file used to generate the annotations! We recommend merging new annotations in on a regular basis for this reason.

### I have feedback and ideas for you!

Wonderful! This is a very new and complex feature; we would LOVE to [hear your feedback](contact) :)
