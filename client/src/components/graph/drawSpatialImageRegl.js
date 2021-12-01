export default function drawSpatialImageRegl(regl) {
  return regl({
    frag: `
        precision mediump float;
        uniform sampler2D texture;
        varying vec2 uv;
        void main () {
          gl_FragColor = texture2D(texture, uv);
        }`,

    vert: `
        precision mediump float;
        attribute vec2 position;
        varying vec2 uv;
        const float zMiddle = 0.;
        uniform mat3 projView;
        
        vec2 norm(vec2 position) {
          return ((position - 0.5) * 2.0) * vec2(1., -1.);
        }
        void main() {
          uv = position;
          vec3 xy = projView * vec3(norm(position), 1.);
          gl_Position = vec4(xy.xy, 0, 1.);
        }`,

    attributes: {
      position: [-1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, 1],
    },

    uniforms: {
      projView: regl.prop("projView"),
      texture: regl.prop("spatialImageAsTexture"),
    },

    count: 6,
  });
}
