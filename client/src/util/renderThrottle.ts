export default function renderThrottle(
  callback: () => void
): (this: unknown) => void {
  /*
  This wraps a call to requestAnimationFrame(), enforcing a single
  render callback at any given time (ie, you can call this any number
  of times, and it will coallesce multiple inter-frame calls into a
  single render).
  */
  let rafCurrentlyInProgress: number | null = null;
  return function f(this: unknown) {
    if (rafCurrentlyInProgress) return;
    // eslint-disable-next-line @typescript-eslint/no-this-alias --- required for functionality
    const context = this;
    rafCurrentlyInProgress = window.requestAnimationFrame(() => {
      callback.apply(context);
      rafCurrentlyInProgress = null;
    });
  };
}
