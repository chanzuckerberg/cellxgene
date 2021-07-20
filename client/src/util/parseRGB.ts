import scaleRGB from "./scaleRGB";

// maintain a cache of already parsed RGB names, as it is reasonably expensive
// to do this operation.  This lets us have speed, but keep the pleasant ability
// to talk about colors by their text description eg, 'rgb(0,0,1)'
//
const colorCache = {};

function parseColorName(c: any) {
  if (c[0] !== "#") {
    const _c = c.replace(/[^\d,.]/g, "").split(",");
    return [scaleRGB(+_c[0]), scaleRGB(+_c[1]), scaleRGB(+_c[2])];
  }
  const parsedHex = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(c);
  return [
    // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
    scaleRGB(parseInt(parsedHex[1], 16)),
    // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
    scaleRGB(parseInt(parsedHex[2], 16)),
    // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
    scaleRGB(parseInt(parsedHex[3], 16)),
  ];
}

export default (c: any) => {
  // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
  let cv = colorCache[c];
  if (!cv) {
    cv = parseColorName(c);
    // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
    colorCache[c] = cv;
  }
  return cv;
};
