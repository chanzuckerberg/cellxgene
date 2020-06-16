export const jestEnv = process.env.JEST_ENV;
export const appPort = process.env.CXG_SERVER_PORT;
export const appUrlBase =
	process.env.CXG_URL_BASE || `http://localhost:${appPort}`;
export const DEV = jestEnv === "dev";
export const DEBUG = jestEnv === "debug";
export const DATASET = "pbmc3k";

if (DEBUG) jest.setTimeout(2 * 60 * 1000);
if (DEV) jest.setTimeout(30 * 1000);
if (!DEBUG && !DEV) jest.setTimeout(10 * 1000);
