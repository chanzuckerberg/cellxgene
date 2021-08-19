import { vec2, mat3 } from "gl-matrix";

const EPSILON = 0.000001;

const scaleSpeed = 0.5;
const scaleMax = 3.0;
const scaleMin = 0.5;
const panBound = 0.8;

// private
const scratch0 = new Float32Array(16);
const scratch1 = new Float32Array(16);

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function clamp(val: any, rng: any) {
  return Math.max(Math.min(val, rng[1]), rng[0]);
}

class Camera {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  canvas: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  prevEvent: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  viewMatrix: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  viewMatrixInv: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  constructor(canvas: any) {
    this.prevEvent = {
      clientX: 0,
      clientY: 0,
      type: 0,
    };
    this.canvas = canvas;
    this.viewMatrix = mat3.create();
    this.viewMatrixInv = mat3.create();
  }

  view() {
    return this.viewMatrix;
  }

  invView() {
    return this.viewMatrixInv;
  }

  distance() {
    return this.viewMatrix[0];
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  pan(dx: any, dy: any) {
    const m = this.viewMatrix;
    const dyRange = [
      -panBound - (m[7] + 1) / m[4],
      panBound - (m[7] - 1) / m[4],
    ];
    const dxRange = [
      -panBound - (m[6] + 1) / m[0],
      panBound - (m[6] - 1) / m[0],
    ];

    const dxClamped = clamp(dx, dxRange);
    const dyClamped = clamp(dy, dyRange);
    if (Math.abs(dxClamped) <= EPSILON && Math.abs(dyClamped) <= EPSILON)
      return;

    mat3.translate(m, m, [dxClamped, dyClamped]);
    mat3.invert(this.viewMatrixInv, m);
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  zoomAt(d: any, x = 0, y = 0) {
    /*
    Camera zoom at [x,y]
    */
    const m = this.viewMatrix;
    const bounds = [-panBound, panBound];
    x = clamp(x, bounds);
    y = clamp(y, bounds);

    const dClamped = clamp(d * m[0], [scaleMin, scaleMax]) / m[0];
    if (Math.abs(1 - dClamped) <= EPSILON) return; // noop request

    mat3.translate(m, m, [x, y]);
    mat3.scale(m, m, [dClamped, dClamped]);
    mat3.translate(m, m, [-x, -y]);

    mat3.invert(this.viewMatrixInv, m);
  }

  /*
  Event handling
  */

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  flush(e: any) {
    this.prevEvent.type = e.type;
    this.prevEvent.clientX = e.clientX;
    this.prevEvent.clientY = e.clientY;
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  localPosition(target: any, canvasX: any, canvasY: any, projectionInvTF: any) {
    /*
    Convert mouse position to local
    */
    const { height, width } = target;
    const targetRect = target.getBoundingClientRect();
    canvasX -= targetRect.left;
    canvasY -= targetRect.top;

    const pos = vec2.fromValues(
      2 * (canvasX / width) - 1,
      -2 * (canvasY / height) + 1
    );
    if (projectionInvTF) {
      vec2.transformMat3(pos, pos, projectionInvTF);
    }
    vec2.transformMat3(pos, pos, this.invView());
    return pos;
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  mousePan(e: any, projectionTF: any) {
    const projectionInvTF = mat3.invert(scratch0, projectionTF);
    const pos = this.localPosition(
      this.canvas,
      e.clientX,
      e.clientY,
      projectionInvTF
    );
    const prev = this.localPosition(
      this.canvas,
      this.prevEvent.clientX,
      this.prevEvent.clientY,
      projectionInvTF
    );

    const delta = vec2.sub(scratch1, pos, prev);
    this.pan(delta[0], delta[1]);
    return true;
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  wheelZoom(e: any, projectionTF: any) {
    const { height } = this.canvas;
    const { deltaY, deltaMode, clientX, clientY } = e;
    const scale = scaleSpeed * (deltaMode === 1 ? 12 : 1) * (deltaY || 0);

    const projectionInvTF = mat3.invert(scratch0, projectionTF);
    const pos = this.localPosition(
      this.canvas,
      clientX,
      clientY,
      projectionInvTF
    );
    this.zoomAt(1 / Math.exp(scale / height), pos[0], pos[1]);
    return true;
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  handleEvent(e: any, projectionTF: any) {
    /*
    process the event, and return true if camera view changed
    */
    let viewChanged = false;
    switch (e.type) {
      case "mousemove": {
        /* eslint-disable-next-line no-bitwise --- MouseEvent.buttons exposes a bitmask, best acted on with bitops */
        if (e.buttons & 0x1) {
          viewChanged = this.mousePan(e, projectionTF);
        }
        this.flush(e);
        break;
      }

      case "wheel": {
        viewChanged = this.wheelZoom(e, projectionTF);
        this.flush(e);
        break;
      }

      default:
        // noop
        break;
    }
    return viewChanged;
  }
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
function attachCamera(canvas: any) {
  return new Camera(canvas);
}

export default attachCamera;
