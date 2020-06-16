export const jestEnv = process.env.JEST_ENV;
export const appPort = process.env.CXG_SERVER_PORT;
export const appUrlBase =
  process.env.CXG_URL_BASE || `http://localhost:${appPort}`;
export const DEV = "dev";
export const DEBUG = "debug";
export const DATASET = "pbmc3k";

export const isDev = jestEnv === DEV;
export const isDebug = jestEnv === DEBUG;
