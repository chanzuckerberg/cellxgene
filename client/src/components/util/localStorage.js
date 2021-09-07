export const KEYS = {
  COOKIE_DECISION: "cxg.cookieDecision",
};

export function storageGet(key, defaultValue = null) {
  try {
    const val = window.localStorage.getItem(key);
    if (val === null) return defaultValue;
    return val;
  } catch (e) {
    return defaultValue;
  }
}

export function storageSet(key, value) {
  try {
    window.localStorage.setItem(key, value);
  } catch {
    // continue
  }
}
