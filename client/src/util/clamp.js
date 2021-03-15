export default function clamp(val, rng) {
  return Math.max(Math.min(val, rng[1]), rng[0]);
}
