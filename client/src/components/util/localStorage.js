export const KEYS = {
  COOKIE_DECISION: "cxg.cookieDecision",
  LOGIN_PROMPT: "cxg.LOGIN_PROMPT",
  WORK_IN_PROGRESS_WARN: "cxg.WORK_IN_PROGRESS_WARN",
};

// TODO(cc) review location
export const WORK_IN_PROGRESS_WARN_STATE = {
  OFF: "off",
  ON: "on",
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
