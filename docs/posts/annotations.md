---
layout: default
title: annotations
description: Creating annotations
---

# Creating annotations in cellxgene

We are _piloting_ a new feature in cellxgene that enables users to create and edit categorical annotations within the app. We'd love for you to try it out and [give us feedback](contact)!

## Quick start for annotations  

You can enable this experimental feature like so:

`cellxgene launch mydata.h5ad --experimental-annotations`

To preserve data provenance, **`cellxgene` does not alter the input h5ad file**.  
Rather, newly-created annotations are saved in a specified CSV file:  
- You will be prompted to enter a name for your annotations the first time you create a new category.  
- We also assign a unique identifier in the form of an 8-character suffix, `########`; this helps cellxgene identify your file to avoid overwriting your work.
- Any annotations you create in the application will be autosaved in `cwd/name-########.csv`, where `cwd` is your current working directory (i.e., the directory you were in when you started cellxgene).

## Data management

### Loading, editing and updating existing draft annotations  

If you'd like to specify the complete file path for your annotations, you can do so by running:  
```
cellxgene launch mydata.h5ad --experimental-annotations-file path/to/myfile.csv
```

If this file already exists and contains compatible annotations, these annotations will be loaded as editable categories that you can update directly. Compatible annotations are tabular, with category names as column headers; `anndata.obs.index` as the index; and categorical values (i.e., fewer unique values per column than specified in `--max-category-items`, default 1000).

Any changes you make will be reflected in the original CSV (which will be overwritten). This is helpful if you wish to annotate over multiple sessions.

If the file does not exist, it will be created.

### Annotations by multiple users  

An alternative to specifying the file path is to specify the output directory, and allow cellxgene to assign filenames. This is most useful for situations where the same cellxgene instance is being used by multiple users to create annotations.

As described in the [hosted](hosted) section, we do not officially support hosted or multi-user use of cellxgene. However, we recognize that the app is often adapted for this purpose, and have tried to provide a "safe path" for multi-user setups that avoids overwriting data.

To specify an output directory, run:  
```
cellxgene launch mydata.h5ad --experimental-annotations-output-dir path/to/annotations-directory/
```

For each user, annotations will be saved as follows:  
- Each user will be prompted to enter a name for their annotations the first time they create a new category.  
- We also assign a unique identifier in the form of an 8-character suffix, `########`; this helps cellxgene identify their specific file to avoid overwriting others' work.
- Any annotations created in the application will be autosaved in `annotations-directory/name-########.csv`


### Merging draft annotations with the main h5ad file

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

### How do you remember my unique ID to match my cellxgene session with my annotations file?
We place a small cookie (file) in your browser that identifies where your draft annotations are saved. This file never leaves your machine, and is never sent to the cellxgene team or anyone else.

### I have feedback and ideas for you!
Wonderful! This is a very new and complex feature; we would _love_ to [hear your feedback](contact) :)
