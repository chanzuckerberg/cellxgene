## Overview

Users may view gene sets along with their dataset and access gene sets created in pip installed cellxgene by following the instructions below.

The order of the gene symbols in each gene set is maintained in the database. 

The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "NOT RECOMMENDED" "MAY", and "OPTIONAL" in this document are to be interpreted as described in BCP 14, [RFC2119](https://www.rfc-editor.org/rfc/rfc2119.txt), and [RFC8174](https://www.rfc-editor.org/rfc/rfc8174.txt) when, and only when, they appear in all capitals, as shown here.

## cellxgene CLI commands

Users MAY use --user-generated-data-dir `file_name` to designate a file where both gene sets and annotations CSVs will be saved.

Users MAY use --gene-sets-file `name_of_file.csv` to designate a CSV with preexisting gene sets you would like to view in cellxgene. The CSV MUST follow the format below.

## cellxgene gene set CSV data format
  
The cellxgene gene set data format is a [*Tidy* CSV](./gene_sets_example.csv) (comma-separated values) file using ASCII encoding. Multiple gene sets MAY be included in the file similar to the [Gene Matrix Transposed](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) or [Gene Matrix](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29) formats.

Example:

| gene_set_name | gene_set_description | gene_symbol | gene_description | provenance1      | provenance1_description |
|---------------|----------------------|-------------|------------------|------------------|-------------------------|
| club.cell     | description          | CCKAR       | description      | Pubmed ID XYZ123 | Primary Pubmed ID       |
| club.cell     |                      | SCGB3A2     | description      | Pubmed ID ABC456 | Primary Pubmed ID       |
| club.cell     |                      | CYP2F2      | description      | Pubmed ID DCF678 | Primary Pubmed ID       |
| macrophage    | description          | CD68        | description      |                  |                         |
| macrophage    |                      | CD163       | description      |                  |                         |


## Mandatory Header

The first row MUST contain header columns using reserved names in the following order:

* `gene_set_name`
* `gene_set_description`
* `gene_symbol`
* `gene_description`

Users MAY include additional header columns. Once gene sets are editted from cellxgene, the additional columns will no longer be stored in the file designated using the command --gene-sets-file.

## Rows

Each subsequent row describes a gene in a gene set.

---

### `gene_set_name`

The `gene_set_name` column MUST contain a value. The values for `gene_set_name` MUST NOT contain illegal ASCII characters or sequences. If the following cases are detected, validation MUST display an error message and fail the upload:

* control characters (decimal 0-31)
* DEL (decimal 127)
* leading spaces (decimal 32) in a field - "     This is an example"
* trailing spaces (decimal 32) in a field - "This is an example     " 
* multiple spaces (decimal 32) "internal" to a field - "This     is an example"

If the `gene_set_name` is missing, validation will display an error message and fail the upload. This is illustrated by **~~?~~** in the example:

| gene_set_name | gene_set_description | gene_symbol | gene_description |
|---------------|----------------------|-------------|------------------|
| club.cell     | description          | CCKAR       | description      |
|     **~~?~~** |                      | SCGB3A2     | description      |
| club.cell     |                      | CYP2F2      | description      |


`gene_symbol(s)` for a `gene_set_name` MAY exist on noncontiguous rows. An out-of-order `gene_symbol` MUST be added to an existing `gene_set_name`. This is illustrated in the example, where **~~CD163~~** is added to the **club.cell** gene set:

| gene_set_name | gene_set_description | gene_symbol | gene_description |
|---------------|----------------------|-------------|------------------|
| club.cell     | description          | CCKAR       | description      |
| club.cell     |                      | SCGB3A2     | description      |
| club.cell     |                      | CYP2F2      | description      |
| macrophage    | description          | CD68        | description      |
| club.cell     |                      | **~~CD163~~**       | description      |

---

### `gene_set_description`

`gene_set_description` MAY contain a value. The first instance of a `gene_set_description` column for a specific `gene_set_name` will be surfaced when a user hovers over `gene_set_name` in cellxgene. All other instances are ignored in subsequent rows for the same `gene_set_name`.

---

### `gene_symbol`

The `gene_symbol` column MUST contain a value that is unique for the `gene_set_name`. Upload will fail and display an error message for a duplicate `gene_symbol`. This is illustrated by **~~CCKAR~~** in the example:

| gene_set_name | gene_set_description | gene_symbol | gene_description |
|---------------|----------------------|-------------|------------------|
| club.cell     | description          | CCKAR       | description      |
| club.cell     |                      | SCGB3A2     | description      |
| club.cell     |                      | **~~CCKAR~~**   | description      |

---

### `gene_description`

The `gene_description` column MAY contain a value and will be surfaced when a user hovers on `gene_symbol` in cellxgene.

## Creating and editting gene sets using cellxgene UX

Users will be able to create and edit gene sets from cellxgene directly. All changes will be saved in the file designated with the command --gene-sets-file. If no file is specified, a CSV file will be created in the user's default directory.

