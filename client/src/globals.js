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
// export const logoColor = "#E9429A";
export const logoColor = "black";

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

export const ordinalColors = [
  "#0ac115",
  "#c10ab6",
  "#c1710a",
  "#0a5ac1",
  "#c1150a",
  "#0ab6c1",
  "#5ac10a",
  "#710ac1",
  "#0ac171",
  "#c10a5a",
  "#b6c10a",
  "#150ac1",
  "#b2ffb7",
  "#ffb2fa",
  "#ffddb2",
  "#b2d4ff",
  "#ffb7b2",
  "#b2faff",
  "#d4ffb2",
  "#ddb2ff",
  "#b2ffdd",
  "#ffb2d4",
  "#faffb2",
  "#b7b2ff",
  "#27a908",
  "#8b08a9",
  "#a93a08",
  "#0877a9",
  "#a90827",
  "#08a98b",
  "#77a908",
  "#3a08a9",
  "#08a93a",
  "#a90877",
  "#a98b08",
  "#0827a9",
  "#00ff0f",
  "#ff00ef",
  "#ff8e00",
  "#0070ff",
  "#ff0f00",
  "#00efff",
  "#70ff00",
  "#8e00ff",
  "#00ff8e",
  "#ff0070",
  "#efff00",
  "#0f00ff",
  "#006606",
  "#66005f",
  "#663900",
  "#002c66",
  "#660600",
  "#005f66",
  "#2c6600",
  "#390066",
  "#006639",
  "#66002c",
  "#5f6600",
  "#060066",
  "#83ff65",
  "#e165ff",
  "#ff9565",
  "#65cfff",
  "#ff6583",
  "#65ffe1",
  "#cfff65",
  "#9565ff",
  "#65ff95",
  "#ff65cf",
  "#ffe165",
  "#6583ff",
  "#009909",
  "#99008f",
  "#995500",
  "#004399",
  "#990900",
  "#008f99",
  "#439900",
  "#550099",
  "#009955",
  "#990043",
  "#8f9900",
  "#090099",
  "#d9fecc",
  "#f1ccfe",
  "#fed7cc",
  "#ccf3fe",
  "#feccd9",
  "#ccfef1",
  "#f3fecc",
  "#d7ccfe",
  "#ccfed7",
  "#feccf3",
  "#fef1cc",
  "#ccd9fe",
  "#47ea51",
  "#ea47e0",
  "#eaa247",
  "#478fea",
  "#ea5147",
  "#47e0ea",
  "#8fea47",
  "#a247ea",
  "#47eaa2",
  "#ea478f",
  "#e0ea47",
  "#5147ea"
];
