import { vec2, mat3 } from "gl-matrix";

const scaleSpeed = 0.5;
const scaleMax = 3.0;
const scaleMin = 0.5;

// private
const scratch0 = new Float32Array(16);
const scratch1 = new Float32Array(16);

class Camera {
  constructor(canvas) {
    this.prevEvent = {
      clientX: 0,
      clientY: 0,
      type: 0
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

  pan(dx, dy) {
    mat3.translate(this.viewMatrix, this.viewMatrix, [dx, dy]);
    mat3.invert(this.viewMatrixInv, this.viewMatrix);
  }

  zoomAt(d, x = 0, y = 0) {
    /*
    Set camera zoom at [x,y]
    */
    const dClamped =
      Math.max(Math.min(d * this.viewMatrix[0], scaleMax), scaleMin) /
      this.viewMatrix[0];
    if (dClamped === 1) return;

    mat3.translate(this.viewMatrix, this.viewMatrix, [x, y]);
    mat3.scale(this.viewMatrix, this.viewMatrix, [dClamped, dClamped]);
    mat3.translate(this.viewMatrix, this.viewMatrix, [-x, -y]);

    mat3.invert(this.viewMatrixInv, this.viewMatrix);
  }

  /*
  Event handling
  */

  flush(e) {
    this.prevEvent.type = e.type;
    this.prevEvent.clientX = e.clientX;
    this.prevEvent.clientY = e.clientY;
  }

  localPosition(target, canvasX, canvasY, projectionInvTF) {
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

  mousePan(e, projectionTF) {
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

  wheelZoom(e, projectionTF) {
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

  handleEvent(e, projectionTF) {
    /*
    process the event, and return true if camera view changed
    */
    let viewChanged = false;
    switch (e.type) {
      case "mousemove": {
        /* eslint-disable no-bitwise */
        if (e.buttons & 0x1) {
          viewChanged = this.mousePan(e, projectionTF);
        }
        /* eslint-enable no-bitwise */
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

function attachCamera(canvas) {
  return new Camera(canvas);
}

export default attachCamera;
