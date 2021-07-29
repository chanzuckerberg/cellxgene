export const KEYS = {
  COOKIE_DECISION: "cxg.cookieDecision",
  LOGIN_PROMPT: "cxg.LOGIN_PROMPT",
};

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function storageGet(key: any, defaultValue = null) {
  try {
    const val = window.localStorage.getItem(key);
    if (val === null) return defaultValue;
    return val;
  } catch (e) {
    return defaultValue;
  }
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function storageSet(key: any, value: any) {
  try {
    window.localStorage.setItem(key, value);
  } catch {
    // continue
  }
}
