import scaleRGB from "./scaleRGB";

// maintain a cache of already parsed RGB names, as it is reasonably expensive
// to do this operation.  This lets us have speed, but keep the pleasant ability
// to talk about colors by their text description eg, 'rgb(0,0,1)'
//
const colorCache = {} as { [key: string]: [number, number, number] };

function parseColorName(c: string): [number, number, number] {
  if (c[0] !== "#") {
    const _c = c.replace(/[^\d,.]/g, "").split(",");
    return [scaleRGB(+_c[0]), scaleRGB(+_c[1]), scaleRGB(+_c[2])];
  }
  const parsedHex = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(c);
  if (!parsedHex || parsedHex.length < 4)
    throw new Error(`Invalid hex color: ${c}`);
  return [
    scaleRGB(parseInt(parsedHex[1], 16)),
    scaleRGB(parseInt(parsedHex[2], 16)),
    scaleRGB(parseInt(parsedHex[3], 16)),
  ];
}

export default (c: string): [number, number, number] => {
  let cv = colorCache[c];
  if (!cv) {
    cv = parseColorName(c);
    colorCache[c] = cv;
  }
  return cv;
};
