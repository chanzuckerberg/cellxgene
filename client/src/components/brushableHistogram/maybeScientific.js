import significantDigits from "../../util/significantDigits";

export default function maybeScientific(x) {
  let format = ",";
  const _ticks = x.ticks(4);

  if (x.domain().some((n) => Math.abs(n) >= 10000)) {
    /* 
          heuristic: if the last tick d3 wants to render has one significant
          digit ie., 2000, render 2e+3, but if it's anything else ie., 42000000 render
          4.20e+n
        */
    format = significantDigits(_ticks[_ticks.length - 1]) === 1 ? ".0e" : ".2e";
  }

  return format;
}
