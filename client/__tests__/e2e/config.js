import * as ENV_DEFAULT from "../../../environment.default.json";

export const jestEnv = process.env.JEST_ENV || ENV_DEFAULT.JEST_ENV;
export const appPort = process.env.APP_PORT || ENV_DEFAULT.CXG_CLIENT_PORT;
export const appUrlBase =
  process.env.CXG_URL_BASE || `http://localhost:${appPort}`;
export const DATASET = "pbmc3k";
export const isDev = jestEnv === ENV_DEFAULT.DEV;
export const isDebug = jestEnv === ENV_DEFAULT.DEBUG;
