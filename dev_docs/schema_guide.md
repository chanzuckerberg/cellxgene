# Cellxgene Schema Guide

Datasets included in the [data portal](https://cellxgene.cziscience.com/) and hosted cellxgene need to follow the schema
described [here](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md). That
schema defines some required fields, requirements about feature labels, and some optional fields that mostly help with
presentation.

The number of fields is rather low, and we expect that information needed to populate those fields should either already
be present in datasets prepared by a submitter or be easy to obtain. However, this still leaves the task of actually
manipulating the dataset so that it follows the schema: adjusting field names, ensuring proper ontologies are used,
converting gene symbols to a common set, etc. This can be tedious and error-prone, and at the beginning of the hosted
cellxgene project, this was always done with engineering support. As we increase the rate at which we add data, we want
to eliminate the need for engineering support so that ultimately submitters themselves can create files that follow the
schema.

## `cellxgene schema apply`

To enable this, we have a new cellxgene subcommand, `cellxgene schema`, that handles applying and verify the schema. Its
first subcommand, `cellxgene schema apply`, takes three inputs:

1. A source h5ad file. The input needs to be an AnnData file, so if a submitter has, say, a serialized Seurat or
   SingleCellExperiment object, it needs to be converted to AnnData first. This can be done with
   (sceasy)[https://github.com/cellgeni/sceasy] or via
   (Seurat)[https://satijalab.org/seurat/v3.1/conversion_vignette.html].
2. A configuration yaml file that describes the fields to add and conversions to apply (see below).
3. A name for the new h5ad file that should follow the schema.

### Configuration yaml

The configuration yaml file describes how to apply the schema. This is an example of a "skeleton" yaml that has all the
fields required for the 1.0.0 schema but is not yet filled in with any logic.

```
uns:
    version:
        corpora_schema_version: 1.0.0
        corpora_encoding_version: 0.1.0
    contributors:
    title:
    layer_descriptions:
    preprint_doi:
    publication_doi:
    organism_ontology_term_id:
obs:
    tissue_ontology_term_id:
    assay_ontology_term_id:
    disease_ontology_term_id:
    cell_type_ontology_term_id:
    sex:
    ethnicity_ontology_term_id:
    development_stage_ontology_term_id:
fixup_gene_symbols:
```

#### Unstructured metadata
The first section is `uns`, which includes metadata fields that describe the whole dataset (see 
(here)[https://anndata.readthedocs.io/en/latest/] for further description of `uns` and `obs`.).

The first line is `version`, which is required for most of our tooling to work. The schema version is set at
1.0.0 in the example above, but of course for future versions that should be changed.

Next is `contributors` which describes who is adding the dataset to the portal. If you consult the schema, you see that
contributors is a list where each element can have `name`, `email`, and `institution`. So when filled out, the
`contributors` field should look like this:

```
contributors:
  - name: Mary B. Scientist
    email: mbs@singlecell.edu
    institution: Single-Cell University
  - name: Robert J. Scientist
    email: rjs@usingle.edu
    institution: University of Single Cell
```

`title` is the name of the dataset, and is just a string that gets displayed in the portal and cellxgene to identify the
dataset.

`layer_descriptions` is free text descriptions of the different
(layers)[https://anndata.readthedocs.io/en/latest/anndata.AnnData.layers.html] of the AnnData file. It should look like
this when complete:
```
layer_descriptions:
  X: CPM and logged
  raw.X: raw
```
Note that one of the layers needs to be "raw", that is, the AnnData file must contain raw counts.

The two DOI fields are optional but can be included if the dataset is associated with a publication or preprint. Note
that the DOI should be a full url:
```
publication_doi: https://doi.org/10.1073%2Fpnas.83.15.5372
```

Finally, the `organism_ontology_term_id` field is the species of the donor organism from the NCBITaxon ontology. The
value for _Homo sapiens_ is `NCBITaxon:9606`:
```
organism_ontology_term_id: NCBITaxon:9606
```
Note that the schema also requires a human-readable `organism` field, but this doesn't need to be included in the yaml.
When the `cellxgene schema apply` script encounters an ontology field, it looks up the label for the term(s) and inserts it
into the appropriate field.


#### Observation metadata
The next section is `obs`, which is metadata than can very for each observation (and "observation" usually means cell).
These fields are all ontology fields except for `sex`, which has its own enumerated set of permitted values.

There are two ways to fill in the `obs` fields. The first is useful when there is only one value for all the
observations in the dataset. This is not uncommon, for example all cells often come from the same assay. In that case
just insert the ontology term:
```
assay_ontology_term_id: EFO:0009922
```

The second is for when there is an existing field in the dataset that needs to be mapped to the schema field. For
example, the submitter may have included cell type annotations in a field called `CellType`, and those annotations may
just be free text. This doesn't follow the schema because it needs to be in `cell_type_ontology_term_id` and
`cell_type`, and it needs ontology terms and labels, not just any text. In that case the field can be a dictionary:

```
cell_type_ontology_term_id:
  CellType:
    t-cell: CL:0000084
    b-cell: CL:0000236
```

This will look at the `obs.CellType` field in the dataset, and where it has the value "T cell", it will insert
`CL:0000084` into `cell_type_ontology_term_id` and its label `T cell` into `cell_type`.

Now there are often situations where there is no valid ontology term for some field. For example, the dataset may have
been produced via an assay not present in `EFO`. Or, a particular cell type may have no entry in `CL`. In that case, a
free text description can be used in the `ontology_term_id` field:

```
assay_ontology_term_id: Sci-Plex
cell_type_ontology_term_id:
  CellType:
    t-cell: CL:0000084
    b-cell: CL:0000236
    new cell type: new cell type
```

In these cases, the `cellxgene schema apply` script will leave the ontology field blank and move the free text
description into the label field. So the `assay_ontology_term_id` in the new dataset would be `""` but `assay` would be
`Sci-Plex`.


#### Gene symbol harmonization

The last section describes how gene symbol conversion should be applied to each of the layers. This is similar to the
`layer_descriptions` field above, but there are only three permitted values: `raw`, `log1p`, and `sqrt`:

```
fixup_gene_symbols:
  X: log1p
  raw.X: raw
```

This tells the script how each each layer was transformed from raw values that can be directly summed. `raw` means that
the layer contains raw counts or some linear tranformation of raw counts. `log1p` means that the layer has `log(X + 1)`
for each the raw `X` values. `sqrt` means `sqrt(X)` (this is not common). For layers produced by Seurat's normalization
or SCTransform functions, the correct choice is usually `log1p`.


### `cellxgene schema validate`

The next `cellxgene schema` subcommand validates that a given h5ad follows a version of the schema. It accepts two
parameters:

1. The h5ad file to check
2. The version of the schema to check against.

If the validation succeeds, the command will have a zero exit code. If it does not, it will have a non-zero exit code
and will print validation failure messages.
