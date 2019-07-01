import { Colors } from "@blueprintjs/core";

/* if a categorical metadata field has more options than this, truncate */
export const maxCategoricalOptionsToDisplay = 100;

/*
these are default values for configuration  the CLI may supply.
See the REST API and CLI specs for more info.
*/
export const configDefaults = {
  features: {},
  displayNames: {},
  parameters: {
    "max-category-items": 1000
  }
};

/* colors */
export const blue = Colors.BLUE3;
export const linkBlue = Colors.BLUE5;
export const lightestGrey = "rgb(249,249,249)";
export const lighterGrey = "rgb(245,245,245)";
export const lightGrey = "rgb(211,211,211)";
export const mediumGrey = "rgb(153,153,153)";
export const darkGrey = "rgb(102,102,102)";
export const darkerGrey = "rgb(51,51,51)";

export const brightBlue = "#4a90e2";
export const brightGreen = "#A2D729";
export const darkGreen = "#448C4D";

export const nonFiniteCellColor = lightGrey;
export const defaultCellColor = "rgb(0,0,0,1)";
export const logoColor = "black"; /* logo pink: "#E9429A" */

/* typography constants */

export const tiniestFontSize = 12;
export const largestFontSize = 24;
export const uppercaseLetterSpacing = "0.04em";
export const bolder = 700;
export const accentFont = "Georgia,Times,Times New Roman,serif";
export const maxParagraphWidth = 600;

/* layout styles constants */

export const cellxgeneTitleLeftPadding = 14;
export const cellxgeneTitleTopPadding = 7;

export const maxControlsWidth = 800;

export const graphMargin = { top: 20, right: 10, bottom: 30, left: 40 };
export const graphWidth = 700;
export const graphHeight = 700;
export const scatterplotMarginLeft = 25;

export const leftSidebarWidth = 365;
export const leftSidebarSectionHeading = {
  fontSize: 18,
  textTransform: "uppercase",
  fontWeight: 500,
  letterSpacing: ".05em"
};
export const leftSidebarSectionPadding = 10;

let _API = {
  // prefix: "http://api.clustering.czi.technology/api/",
  // prefix: "http://tabulamuris.cxg.czi.technology/api/",
  // prefix: "http://api-staging.clustering.czi.technology/api/",
  prefix: "http://localhost:5005/api/",
  version: "v0.2/"
};

if (window.CELLXGENE && window.CELLXGENE.API) _API = window.CELLXGENE.API;
export const API = _API;
