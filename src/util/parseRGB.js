import {scaleRGB} from "./scaleRGB";

export const parseRGB = (c) => {
  if (c[0] !== "#") {
    const _c = c.replace(/[^\d,.]/g, '').split(',');
    return [
      scaleRGB(+_c[0]),
      scaleRGB(+_c[1]),
      scaleRGB(+_c[2])
    ];
  } else {
    var parsedHex = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(c);
    return [
      scaleRGB(parseInt(parsedHex[1], 16)),
      scaleRGB(parseInt(parsedHex[2], 16)),
      scaleRGB(parseInt(parsedHex[3], 16))
    ];
  }
};
