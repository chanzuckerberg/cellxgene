export const save = (key, data) => {
  try {
    localStorage.setItem(key, JSON.stringify(data));
  } catch (err) {
    console.warn(`Failed to save ${key} to localStorage: ${err}`);
  }
};

export const load = (key, defaultValue) => {
  try {
    return JSON.parse(localStorage.getItem(key)) || defaultValue;
  } catch (err) {
    console.warn(`Failed to load ${key} from localStorage: ${err}`);
    return defaultValue;
  }
};
