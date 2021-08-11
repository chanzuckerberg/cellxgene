import { Colors } from "@blueprintjs/core";
import { dispatchNetworkErrorMessageToUser } from "./util/actionHelpers";
import ENV_DEFAULT from "../../environment.default.json";

/* overflow category values are created  using this string */
export const overflowCategoryLabel = ": all other labels";

/* default "unassigned" value for user-created categorical metadata */
export const unassignedCategoryLabel = "unassigned";

/*
these are default values for configuration the CLI may supply.
See the REST API and CLI specs for more info.
*/
export const configDefaults = {
  features: {},
  displayNames: {},
  parameters: {
    "disable-diffexp": false,
    "diffexp-may-be-slow": false,
  },
  links: {},
};

/*
Most configuration is stored in the reducer.  A handful of values
are global and stored here.  They are typically set by the config
action handler, which pull the information from the backend/CLI.

All should be set here to their default value.
*/
export const globalConfig = {
  /* if a categorical metadata field has more options than this, truncate */
  maxCategoricalOptionsToDisplay: 200,
};

/* colors */
export const blue = Colors.BLUE3;
export const linkBlue = Colors.BLUE5;
export const lightestGrey = "rgb(249,249,249)";
export const lighterGrey = "rgb(245,245,245)";
export const lightGrey = Colors.LIGHT_GRAY1;
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

export const datasetTitleMaxCharacterCount = 25;

export const maxControlsWidth = 800;

export const graphMargin = { top: 20, right: 10, bottom: 30, left: 40 };
export const graphWidth = 700;
export const graphHeight = 700;
export const scatterplotMarginLeft = 11;

export const rightSidebarWidth = 365;
export const leftSidebarWidth = 365;
export const leftSidebarSectionHeading = {
  fontSize: 18,
  textTransform: "uppercase",
  fontWeight: 500,
  letterSpacing: ".05em",
};
export const leftSidebarSectionPadding = 10;
export const categoryLabelDisplayStringLongLength = 27;
export const categoryLabelDisplayStringShortLength = 11;
export const categoryDisplayStringMaxLength = 33;

export const maxUserDefinedGenes = 25;
export const maxGenes = 100;

export const diffexpPopNamePrefix1 = "Pop1 high";
export const diffexpPopNamePrefix2 = "Pop2 high";

/* various timing-related behaviors */
export const tooltipHoverOpenDelay = 1000; /* ms delay before a tooltip displays */
export const tooltipHoverOpenDelayQuick = 500;

const CXG_SERVER_PORT =
  process.env.CXG_SERVER_PORT || ENV_DEFAULT.CXG_SERVER_PORT;

let _API;

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
if ((window as any).CELLXGENE && (window as any).CELLXGENE.API) {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  _API = (window as any).CELLXGENE.API;
} else {
  if (CXG_SERVER_PORT === undefined) {
    const errorMessage = "Please set the CXG_SERVER_PORT environment variable.";
    dispatchNetworkErrorMessageToUser(errorMessage);
    throw new Error(errorMessage);
  }

  _API = {
    // prefix: "http://api.clustering.czi.technology/api/",
    // prefix: "http://tabulamuris.cxg.czi.technology/api/",
    // prefix: "http://api-staging.clustering.czi.technology/api/",
    prefix: `http://localhost:${CXG_SERVER_PORT}/api/`,
    version: "v0.2/",
  };
}

/*
 TODO(cc) temp set local flag to handle differences between local and deployed environments
 */
_API.local = window.location.hostname === "localhost";

/*
 TODO(cc) temp set Portal staging prefix to handle meta and collection API requests
 */
_API.portalPrefix =
  "https://api.cellxgene.staging.single-cell.czi.technology/dp/v1/";

/*
 TODO(cc) temp set of Portal/Explorer origin, required for breadcrumb links as well as generating explore URL param for
  meta endpoint in environments where hosted origin does not match Portal/dataset deployment URL origin (eg local and
  canary).
 */
_API.origin = "https://cellxgene.staging.single-cell.czi.technology/";

export const API = _API;
