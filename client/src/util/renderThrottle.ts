// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export default function renderThrottle(callback: any) {
  /*
  This wraps a call to requestAnimationFrame(), enforcing a single
  render callback at any given time (ie, you can call this any number
  of times, and it will coallesce multiple inter-frame calls into a
  single render).
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  let rafCurrentlyInProgress: any = null;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  return function f(this: any) {
    if (rafCurrentlyInProgress) return;
    // eslint-disable-next-line @typescript-eslint/no-this-alias --- FIXME: disabled temporarily on migrate to TS.
    const context = this;
    rafCurrentlyInProgress = window.requestAnimationFrame(() => {
      callback.apply(context);
      rafCurrentlyInProgress = null;
    });
  };
}
