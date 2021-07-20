/**
 * The function will ensure a number is above or below given thresholds
 * @param val - a number
 * @param rng - an array of two numbers, a min and a max
 * Ie., in the case of a histogram brush selection:
 * const x0 = x(clamp(selectionRange[0], [min, max]));
 * const x1 = x(clamp(selectionRange[1], [min, max]));
 * @returns a number
 */

export default function clamp(val, rng) {
  return Math.max(Math.min(val, rng[1]), rng[0]);
}
