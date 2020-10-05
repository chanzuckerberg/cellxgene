## UPDATE (9/30/2020): Starting today, the name Corpora will only be used as the internal project name, with cellxgene Data Portal being the official product name

# CXG Data Format Specification

Document Status: _draft_

Version: 0.2.0 (_DRAFT, not yet approved_)

Date Last Modified: 2020-07-23

## Introduction

CXG is a cellxgene-private data format, used for at-rest storage of annotated matrix data. It is similar to [AnnData](https://anndata.readthedocs.io/en/stable/), but with performance and access characteristics amenable to a multi-dataset, multi-user serving environment.

CXG is built upon the [TileDB](https://tiledb.com/) embedded database. Each CXG is a TileDB [group](https://docs.tiledb.com/main/api-usage/object-management), which in turn includes one or more TileDB multi-dimensional arrays.

This document presumes familiarity with [TileDB terminology and concepts](https://docs.tiledb.com/main/), the [Corpora schema](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md) and its [H5AD encoding](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md), and the AnnData/H5AD data model.

This document also leverages the current cellxgene schema, which is documented in the [REST API spec](./REST_API.md).

### Terminology

Unless explicitly noted, the AnnData conventions and terminology are adopted when referring to general annotated matrix characteristics (eg, `n_obs` is the number of observations/rows/cells in the annotated matrix).  Where implied by context, eg, "TileDB array attribute", domain-specific terms are used.

Where capitalized, [IETF RFC 2119](https://www.ietf.org/rfc/rfc2119.txt) conventions are followed (ie, conventions MUST be followed).

_Author's note:_ if you see any ambiguous terms, please call them out for clarification.

### Reserved

The `cxg` prefix is used for CXG-specific names.

### Encoding Data With TileDB Arrays

The TileDB array schema authoritatively defines the characteristics of each array (eg, the type of `X` is defined by [`X.schema`](https://tiledb-inc-tiledb-py.readthedocs-hosted.com/en/stable/python-api.html#tiledb.libtiledb.Array.schema)). In some cases, additional metadata is required for the CXG, and is attched to the array using the TileDB [array metadata](https://docs.tiledb.com/main/basic-concepts/array-metadata) capability.

All TileDB arrays MUST have a uint32 domain, zero based.  All X counts and embedding coordinates SHOULD be coerced to float32, which is ample precision for visualization purposes, and MUST be a numeric type. Dataframe (metadata) types are generally preserved, or where that is not possible, converted to something with equal representative value in the cellxgene application (eg, categorical types are converted to string, bools to uint8, etc).

CXG consumers (readers) MUST be prepared to handle any legal TileDB compression, global layout and tile size.  CXG writers SHOULD attempt to encode data using best-effort heuristics for time and space considerations (eg, dense/sparse encoding tradeoffs).

## Entities

### CXG

The CXG is a TileDB group containing all data and metadata for a single annotated matrix.  The following objects MUST be present in a CXG, except where noted as optional:
* __obs__: a TileDB array, of shape (n_obs,), containing obs annotations, each annotation stored in a separate TileDB array attribute.
* __var__: a TileDB array, of shape (n_var,), containing var annotations, each annotation stored in a separate TileDB array attribute.
* __X__: a TileDB array, of shape (n_obs, n_var), with a single TileDB attribute of numeric type.
* __X_col_shift__: (optional) TilebDB Array used in column shift encoding, shape (n_var,), dtype = X.dtype. Single unnamed numeric attribute.  
* __emb__: a TileDB group, which in turn contains all (zero or more) embeddings.
* __emb__/\<embedding_name\>__: a TileDB array, with a single anonymous attribute, of numeric type, and shape (n_obs, N>=2).
* __cxg_group_metadata__: an empty TileDB array, used to store CXG-wide metadata

### obs and var

All per-observation (obs) and per-feature (var) data is encoded in a TileDB array named `obs` and `var` respectively, with shape (n_obs,) and (n_var,).  Each TileDB array has an array attribute for each obs/var column.  All TileDB array attributes will have the same type and value as the original data, eg, float32, with the following exceptions:
* bool is encoded as uint8 (1/0)
* categorical is encoded as string
* Numeric types are cast to 32-bit equivalents

In addition to the obs/var data, both TileDB arrays contain an optional 'cxg_schema' metadata field that is a JSON string containing per-column (attribute) schema hinting.  This is used where the TileDB native typing information is insufficient to reconstruct useful information such as categorical typing from Pandas DataFrames, and to communicate which column is the preferred human-readable index for obs & var.

The `cxg_schema` JSON string is attached to the TileDB array metadata, and is a dictionary containing the following top-level names:
* "index": string, containing the name of the index column
* \<column-name\>: optional, a JSON dict, contain a schema definition using the same format as the cellxgene REST API /schema route

For example:
```
{
    "index": "obs_index",
    "louvain": { "type": "categorical", "categories": [ "0", "1", "2", "3", "4" ]}
    "is_useful": { "type": "boolean" }
}
```

### X

TileDB array, with a single anonymous attribute, shape (n_obs, n_var), containing the count matrix (equivalent to the AnnData `X` array).  MUST have numeric type, and SHOULD be float32.  The TileDB schema defines type and sparsity, and both dense and sparse encoding are supported.

### X_col_shift

Optional TileDB array, used to encode-per column offsets for column-shift sparse encoding.  The TileDB array will have a single anonymous attribute, of the same type as the X array, and shape (n_var,).  

If the X array is sparse, and X_col_shift exists, then all values in the i'th column were subtracted by X_col_shift[i].

### emb and embedding arrays

A CXG must have a group named `emb`, which will contain all embeddings.  Embeddings are encoded as TileDB arrays, of numeric type and shape (n_obs, >=2).  The arrays SHOULD be coerced to float32, and MUST be a numeric type.  The TileDB array name will be assumed to be the embedding name (conventionally, embedding names in CXG are _not_ prefixed with an `X_` as they are in AnnData).

CXG supports zero or more embeddings. Note that cellxgene currently _requires_ at least one embedding.

### cxg_group_metadata

Required, but empty TileDB array, used to store CXG-wide metadata.  The following fields are defined:
* __cxg_version__: (required) a semver string identifying the specification version used to encode the CXG.
* __cxg_properties__: (optional) a dictionary containing dataset wide properties, defined below.
* __cxg_category_colors__: (optional) a categorical color table, defined below.

#### cxg_properties

The properties metadata dictionary contains dataset-wide properties, encoded as a JSON dictionary. Currently, the following fields are defined:
* title: string, dataset human name (eg, "Lung Tissue")
* about: string, fully-qualified http/https URL, linking to more information on the dataset.

All implementions MUST ignore unrecognized fields. 

#### cxg_category_colors

This optional field contains a copy of the category color table, which MAY be used to display category-specific color labels. This is a JSON dictionary, containing a per-category color-table.  Each color table is named `{category_name}_colors`, and is itself a dictionary mapping label name to RGB color.  For example:

```
{
    "louvain_colors": {
        "0": "#FFFFFF",
        "1": "#000000"
    }
}
```

## Corpora Schema Encoding

The [Corpora schema](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md) and [Corpora AnnData encoding](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md) define a set of metadata and encoding conventions for annotated matrices.  When a Corpora dataset is encoded as a CXG, the following shall apply.

### Corpora metadata property

A CXG containing a Corpora dataset will contain a property in the __cxg_group_metadata__ field named `corpora`. The value will be a JSON encoded string, which in turn contains all properties defined in the [Corpora AnnData uns](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md#uns) container. For example:

```
{
    "corpora": {
        "version": {
            "corpora_schema_version": "1.0.0",
            "corpora_encoding_version": "0.1.0",
        }
    }
}
```

The `corpora` metadata, if present, MUST contain the version information. Optionality of other values in this object will follow the specifications set forth in the relevant Corpora schema specification (ie, optional fields are optional, required are present, etc), with the following changes:
* the contents of `corpora_encoding_version` MUST be identical to the `cxg_version`, as this field is defined as the current object encoding version, *NOT* the source data encoding version.
* the entire encoding will be JSON, rather than a hybrid Python/JSON encoding, but will otherwise follow the data structure defined by the AnnData Corpora encoding.
* the `<obs_column>_colors` will be omitted in favor of `cxg_category_colors`

### Other Corpora fields

All other Corpora schema fields will be encoded into a CXG using the conventions defined in the [Corpora Schema AnnData Implementation](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md). For example, fields in `AnnData.obs` will be encoded in the CXG `obs`array as defined [above](#obs-and-var).

### Compatibility with CXG 0.1.0

For backwards compatibility and continuity with CXG version 0.1.0, the following MUST be implemented.

#### Presentation Hints
* The [Corpora `title`](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md#presentation-metadata) value MUST be saved in the `cxg_properties.title` field.
* The [Corpora `color_map`](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md#presentation-hints), when present in the dataset, MUST be saved in the `cxg_category_colors` field and NOT in the `corpora` field.
* The [Corpora SUMMARY `project_link`](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema.md#presentation-hints), if present, MUST be saved in the `cxg_properties.about` field.

Where these values differ in the final CXG, the `cxg_properties` values WILL take precedence.

## CXG Version History

There were several ad hoc version of CXG created prior to this spec.  This describes the _proposed_ next version of CXG, which incoporates support for Corpora schema semantics.  Prior verisons:
* _unnamed_ - an unnamed development version. Did not include explicit versioning support in the data model, but can be detected by the absence of __cxg_group_metadata__ and any version property. Created in early 2020, and not actively used in production
* 0.1 - the first and current version, defined to support the capabilities of the mid-2020 cellxgene.  Created in early 2020, and in active use. Includes everything in this spec, excluding Corpora schema support.  __NOTE:__ this version is encoded with a short-hand (malformed) semver version number.
* 0.2.0 - this specification.
