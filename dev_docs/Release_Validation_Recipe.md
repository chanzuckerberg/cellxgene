# Pre-release Validation Plan
_This document contains a general purpose test plan for validating cellxgene prior to release_.

All steps are expected to pass with no errors or malfunctions.   Tester should check:
* CLI for errors (eg, build error, engine error)
* Browser console for errors
* Expected CLI and/or browser UI function

## [Backend Test](#backend-test)

### Getting started build & install validation:
Goal: validate build, install and example demo correctness, per Getting Started instructions.

1. Clean clone of the cellxgene repo into a local directory
2. Follow build & install instructions from the Getting Started guide
    * Confirm no build or install errors
    * Confirm `cellxgene --help` functions correctly
3. Follow the example data set demo startup from Getting Started guide and confirm front-end data loads correctly
4. Confirm all package (and other) version dependencies are correct and match 

### Tabula Muris 
Goal: basic functional validation of b/e functions using Tabula Muris data and the scanpy engine.

1. Precondition: cellxgene built & installed.
2. Load the Tabula Muris data set: `cellxgene --title 'T. Muris' scanpy directory-name/`
3. Verify all metadata selectors display and have correct type/options:
    * *TODO:* _need list of metadata and their type_
    * ...
4. Verify default graph display has expected layout.   _need screen shot of expected layout_
5. Verify selection controls work as expected:
    * continuous metadata field
    * categorical metadata field
    * graph/lasso select
6. Verify color by metadata type
7. Select two cell sets and confirm differential expression compute succeeds
6. Verify expression scatter plot is correct 


## [Front-end Compatibility Test](#frontend-compatability-test)
Goal: verify front-end UI compatibility with a given browser variant/version/platform.

1. Start back-end on PBM3K data set
2. Load UI
3. Verify all major UI modes/functions:
    * Graph display
    * Metadata selector display
    * Title display
    * Selection - single and multiple fields - correctly display in cluster graph
    * Regraph & reset function correctly
    * Differential expression calc & scatter plot display
    * All graphs maintain consistent selection state
    * All selection widgets (eg, continuous metadata selector) maintain correct status (consistent with graph displays)
4. Verify overall performance is reasonable/interactive
5. Verify no errors on CLI or browser console


## [End-to-end Functional Test](#e2e-functional-test)
Goal: confirm end-to-end functional behavior is as expected.

### ScanPy engine

1. Basic data load and display
    * ...
2. Select & multi-select of metadata and coordinates
    * ...
3. Color by metadata
    * ...
4. Regraph / reset
    * ...
5. Differential expression: scatterplot, top-N genes, etc.
    * ...
6. Color by expression
    * ...
7. Arbitrary gene expression
    * ...
