---
layout: default
title: annotations
description: Creating annotations
---

# Creating annotations in cellxgene

We are _piloting_ a new feature in cellxgene that enables users to create and edit categorical annotations within the app. We'd love for you to try it out and [give us feedback](contact)!

# Data lifecycle for annotations  
## 1. Creating annotations (quickstart)

To get started, run:
```
cellxgene launch mydata.h5ad
```

To preserve data provenance, **`cellxgene` does not alter the input h5ad file**.  Rather, newly-created annotations are saved in a specified CSV file:  
```
annotations-directory/name-########.csv
```

- The default `annotations-directory` is the directory where `mydata.h5ad` is saved; if we can't find that directory (e.g., you're loading an h5ad from a URL), it falls back to your current working directory (i.e., the directory you were in when you started cellxgene).
- You will be prompted to enter a name for your annotations the first time you create a new category.  
- We also assign a unique identifier in the form of an 8-character suffix, `########`; this helps cellxgene identify your file to avoid overwriting your work.

## 2. Loading, editing and updating existing draft annotations  

Cellxgene allows you to load and edit compatible draft annotations across multiple sessions.

Compatible annotations are tabular, with category names as column headers; `anndata.obs.index` as the index; and categorical values (i.e., fewer unique values per column than specified in `--max-category-items`, default 1000).

There are two options for updating draft annotations.

### Autodetect annotations csv

Cellxgene will automatically find and reload your draft annotations in editable mode.

This assumes that:
1 - The h5ad filename is the same  
2 - You launch cellxgene from the `annotations-directory` (i.e., the directory that contains your CSV)  
3 - You use the same browser and have not cleared your cookies (we use a small cookie to keep track of which user created the file to avoid accidental overwrites; see FAQ)


### Specify an annotations csv  
**This mode is only appropriate for single-user, local cellxgene instances**

If you'd like to specify the complete file path for your annotations, you can do so by running:  
```
cellxgene launch mydata.h5ad --annotations-input-file path/to/myfile.csv
```

Any changes you make will be reflected in the original CSV. If the file does not exist, it will be created.  
**Please note that this file will be overwritten, making this mode inappropriate for hosted / multi-user settings (see below).**

### 3. Merging draft annotations with the main h5ad file

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

# Annotations by multiple users  

As described in the [hosted](hosted) section, we do not officially support hosted or multi-user use of cellxgene. However, we recognize that the app is often adapted for this purpose, and have tried to provide a "safe path" for multi-user setups that avoids overwriting data.

Specifying a single file name for multiple contributors will result in data overwriting. To avoid this, you can instead specify an output directory and allow cellxgene to assign filenames.

To specify an output directory, run:  
```
cellxgene launch mydata.h5ad --annotations-output-dir path/to/annotations-directory/
```

For each user, annotations will be saved as follows:  
- Each user will be prompted to enter a name for their annotations the first time they create a new category.  
- We also assign a unique identifier in the form of an 8-character suffix, `########`; this helps cellxgene identify their specific file to avoid overwriting others' work.
- Any annotations created in the application will be autosaved in `annotations-directory/name-########.csv`



## FAQ

### How do I know my annotations are saved?
`cellxgene` autosaves any changes made to your annotations every 3 seconds.

### I think I deleted my annotations! Oh noes!  
Not to worry! We save the last 10 versions of your annotations in `annotations-directory/NAME-backups/`

### What about creating continuous annotations?  
Continuous metadata is important! However, these values (e.g., pseudotime) are the result of statistical analyses that are beyond cellxgene's visualization- and exploration-focused scope. We do, of course, provide visualization of continuous metadata values computed elsewhere and stored in `anndata.obs`.

### I keep getting weird index errors when trying to join my annotations to my anndata??  
This is most likely because the h5ad file you are working with is not the original file used to generate the annotations! We recommend merging new annotations in on a regular basis for this reason.

### How do you remember my unique ID to match my cellxgene session with my annotations file?
We place a small cookie (file) in your browser that identifies where your draft annotations are saved. This file never leaves your machine, and is never sent to the cellxgene team or anyone else.

### I have feedback and ideas for you!
Wonderful! This is a relatively new feature; we would _love_ to [hear your feedback](contact) :)
