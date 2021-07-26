export default function renderThrottle(callback: any) {
  /*
  This wraps a call to requestAnimationFrame(), enforcing a single
  render callback at any given time (ie, you can call this any number
  of times, and it will coallesce multiple inter-frame calls into a
  single render).
  */
  let rafCurrentlyInProgress: any = null;
  return function f(this: any) {
    if (rafCurrentlyInProgress) return;
    const context = this;
    rafCurrentlyInProgress = window.requestAnimationFrame(() => {
      callback.apply(context);
      rafCurrentlyInProgress = null;
    });
  };
}
