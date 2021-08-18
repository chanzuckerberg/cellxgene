export default function renderThrottle<T>(
  callback: (this: T) => void
): (this: T) => void {
  /*
  This wraps a call to requestAnimationFrame(), enforcing a single
  render callback at any given time (ie, you can call this any number
  of times, and it will coallesce multiple inter-frame calls into a
  single render).
  */
  let rafCurrentlyInProgress: number | null = null;
  return function f(this: T) {
    if (rafCurrentlyInProgress) return; // eslint-disable-next-line @typescript-eslint/no-this-alias --- required for functionality
    rafCurrentlyInProgress = window.requestAnimationFrame(() => {
      callback.call<T, never, void>(this);
      rafCurrentlyInProgress = null;
    });
  };
}
