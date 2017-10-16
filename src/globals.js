/* these will be either (preferably) specified or inferred */
export const categories = ["Sample.type", "Selection", "Location", "Sample.name", "Class", "Neoplastic"];
export const continuous = [
  "Total_reads",
  "Unique_reads",
  "Unique_reads_percent",
  "ERCC_reads",
  "Non_ERCC_reads",
  "ERCC_to_non_ERCC",
  "Genes_detected",
  "Multimapping_reads_percent",
  "Splice_sites_AT.AC",
  "Splice_sites_Annotated",
  "Splice_sites_GC.AG",
  "Splice_sites_GT.AG",
  "Splice_sites_non_canonical",
  "Splice_sites_total",
  "Unmapped_mismatch",
  "Unmapped_other",
  "Unmapped_short"
]

/* colors */
export const hcaBlue = "#1c7cc7"
export const lighterGrey = "rgb(245,245,245)";
export const lightGrey = "rgb(211,211,211)";
export const mediumGrey = "rgb(153,153,153)";
export const darkGrey = "rgb(102,102,102)";
export const darkerGrey = "rgb(51,51,51)";

export const tiniestFontSize = 12;


export const bolder = 700;

export const API = {
  prefix: "http://api.clustering.czi.technology/api/",
  version: "v0.1/",
}

export const accentFont = "Georgia,Times,Times New Roman,serif";
export const maxParagraphWidth = 600;
export const maxControlsWidth = 800;
