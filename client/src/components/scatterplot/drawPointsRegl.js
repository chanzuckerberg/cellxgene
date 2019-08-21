// jshint esversion: 6
const mat4 = require("gl-mat4");

export default function(regl) {
  return regl({
    vert: `
    precision mediump float;
    attribute vec2 position;
    attribute vec3 color;
    attribute float size;
    uniform mat3 projection;
    varying vec3 fragColor;
    void main() {
      gl_PointSize = 0.5 + size;
      vec3 xy = projection * vec3(position, 1.);
      gl_Position = vec4(xy.xy, 0., 1.);
      fragColor = color;
    }`,

    frag: `
    precision mediump float;
    varying vec3 fragColor;
    void main() {
      if (length(gl_PointCoord.xy - 0.5) > 0.5) {
        discard;
      }
      gl_FragColor = vec4(fragColor, 1.);
    }`,

    attributes: {
      position: regl.prop("position"),
      color: regl.prop("color"),
      size: regl.prop("size")
    },

    uniforms: {
      projection: regl.prop("projection")
    },

    count: regl.prop("count"),

    primitive: "points"
  });
}
