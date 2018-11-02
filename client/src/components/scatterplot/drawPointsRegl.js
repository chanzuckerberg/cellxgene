// jshint esversion: 6
const mat4 = require("gl-mat4");

// opacity: https://github.com/spacetx/starfish/blob/master/viz/draw/regions.js

export default function(regl) {
  return regl({
    vert: `
    precision mediump float;
    attribute vec2 position;
    attribute vec3 color;
    attribute float size;
    uniform float distance;
    uniform mat4 projection, view;
    varying vec3 fragColor;
    void main() {
      gl_PointSize = 7.0 / pow(distance, 2.5) + size;
      gl_Position = projection * view * vec4(position.x, -position.y, 0, 1);
      fragColor = color;
    }`,

    frag: `
    precision mediump float;
    varying vec3 fragColor;
    void main() {
      if (length(gl_PointCoord.xy - 0.5) > 0.5) {
        discard;
      }
      gl_FragColor = vec4(fragColor, 1);
    }`,

    attributes: {
      position: regl.prop("position"),
      color: regl.prop("color"),
      size: regl.prop("size")
    },

    uniforms: {
      distance: regl.prop("distance"),
      view: regl.prop("view"),
      projection: (context, props) =>
        mat4.perspective(
          [],
          Math.PI / 2,
          (context.viewportWidth * props.scale) / context.viewportHeight,
          0.01,
          1000
        )
    },

    count: regl.prop("count"),

    primitive: "points"
  });
}
