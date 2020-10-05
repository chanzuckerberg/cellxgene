# Design Principles

cellxgene is a tool for scientists investigating single cell rna seq datasets. The design of cellxgene proceeds from a small number of core principles:

1. the tool should be highly **scaleable** and handle *exploration* of millions of cells in the browser at interactive speeds 
2. the tool should be **data dense** and give scientists powerful views into data
3. the tool should be powerfully expressive and **optimized for the `nth` day of use** rather than the first day, in the spirit of enterprise tools, even if that requires training or onboarding 
    - (think: photoshop, a broad array of youtube tutorials, books and trainings)
4. the tool should **avoid duplicating data** onscreen 
    - for example, if we render a categorical label on the left hand side bar, locate further information related to that label in place on the left hand sidebar, rather than rendering that data again.
    - for a concrete example, see the relationship between this solution: https://github.com/chanzuckerberg/cellxgene/pull/827 and this problem: https://github.com/chanzuckerberg/cellxgene/issues/762
5. the tool should **enable interactive recomputation of views** into the data, for example: 
    - crossfiltering
    - dataframe subsetting
    - interactive differential expression
    - interactive reprojection of umap
6. the tool should enable a productive workflow between those who are computational and those who are not
