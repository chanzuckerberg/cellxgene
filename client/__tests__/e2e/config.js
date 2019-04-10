export const jest_env = process.env.JEST_ENV || "dev";
export const appPort = process.env.JEST_CXG_PORT || 3000;
export const appUrlBase = `http://localhost:${appPort}`;
export const DEV = jest_env === "dev";
export const DEBUG = jest_env === "debug";
export const DATASET = "pbmc3k";
