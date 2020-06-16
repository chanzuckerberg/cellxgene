export const DEV = "dev";
export const DEBUG = "debug";

export const PROD = "prod";

export const CXG_SERVER_PORT = 5005;

export const CXG_CLIENT_PORT = 3000;

export const jestEnv = process.env.JEST_ENV || PROD;
export const appPort = process.env.CXG_CLIENT_PORT || CXG_CLIENT_PORT;
export const appUrlBase =
  process.env.CXG_URL_BASE || `http://localhost:${appPort}`;

export const DATASET = "pbmc3k";

export const isDev = jestEnv === DEV;
export const isDebug = jestEnv === DEBUG;
