export default function drawSpatialImageRegl(regl) {
  return regl({
    frag: `
        precision mediump float;

        // our texture
        uniform sampler2D u_image;

        // the texCoords passed in from the vertex shader.
        varying vec2 v_texCoord;

        void main() {
            gl_FragColor = texture2D(u_image, v_texCoord);
        }`,

    vert: `
        attribute vec2 a_position;
        attribute vec2 a_texCoord;

        uniform vec2 u_resolution;

        uniform mat3 projView;

        varying vec2 v_texCoord;

        void main() {
            // convert the rectangle from pixels to 0.0 to 1.0
            vec3 pos = vec3(a_position, 1.);
            vec2 zeroToOne = pos.xy / u_resolution;

            // convert from 0->1 to 0->2
            vec2 zeroToTwo = zeroToOne * 2.0;

            // convert from 0->2 to -1->+1 (clipspace)
            vec2 clipSpace = zeroToTwo - 1.0;

            vec3 pos2 = projView * vec3(clipSpace, 1.);

            gl_Position = vec4(pos2.xy , 0, 1);

            // pass the texCoord to the fragment shader
            // The GPU will interpolate this value between points.
            v_texCoord = a_texCoord;
        }`,

    attributes: {
      a_texCoord: [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0],
      // a_position: [
      //     10, 0,
      //     10 + regl.prop("img_width"), 0,
      //     10, 0 + regl.prop("img_height"),
      //     10, 0 + regl.prop("img_height"),
      //     10 + regl.prop("img_width"), 0,
      //     10 + regl.prop("img_width"), 0 + regl.prop("img_height"),
      //     ],
      a_position: [
        0,
        0,
        0 + 1921,
        0,
        0,
        0 + 2000,
        0,
        0 + 2000,
        0 + 1921,
        0,
        0 + 1921,
        0 + 2000,
      ],
    },

    uniforms: {
      projView: regl.prop("projView"),
      u_image: regl.prop("spatialImageAsTexture"),
      color: [1, 0, 0, 1],
      u_resolution: [1921, 2000],
      //   translate:
    },

    count: 6,
  });
}
