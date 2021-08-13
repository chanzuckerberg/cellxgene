import { vec2, mat3 } from "gl-matrix";
import clamp from "./clamp";

const EPSILON = 0.000001;

const scaleSpeed = 0.5;
const scaleMax = 3.0;
const scaleMin = 0.5;
const panBound = 0.8;

// private
const scratch0 = new Float32Array(16);
const scratch1 = new Float32Array(16);

class Camera {
  canvas: HTMLCanvasElement;

  prevEvent: {
    clientX: number;
    clientY: number;
    type: string;
  };

  viewMatrix: mat3;

  viewMatrixInv: mat3;

  constructor(canvas: HTMLCanvasElement) {
    this.prevEvent = {
      clientX: 0,
      clientY: 0,
      type: "",
    };
    this.canvas = canvas;
    this.viewMatrix = mat3.create();
    this.viewMatrixInv = mat3.create();
  }

  view(): mat3 {
    return this.viewMatrix;
  }

  invView(): mat3 {
    return this.viewMatrixInv;
  }

  distance(): number {
    return this.viewMatrix[0];
  }

  pan(dx: number, dy: number): void {
    const m = this.viewMatrix;
    const dyRange: [number, number] = [
      -panBound - (m[7] + 1) / m[4],
      panBound - (m[7] - 1) / m[4],
    ];
    const dxRange: [number, number] = [
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

  zoomAt(d: number, x = 0, y = 0): void {
    /*
    Camera zoom at [x,y]
    */
    const m = this.viewMatrix;
    const bounds: [number, number] = [-panBound, panBound];
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

  flush(e: MouseEvent) {
    this.prevEvent.type = e.type;
    this.prevEvent.clientX = e.clientX;
    this.prevEvent.clientY = e.clientY;
  }

  localPosition(
    target: HTMLCanvasElement,
    canvasX: number,
    canvasY: number,
    projectionInvTF: mat3
  ): vec2 {
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

  mousePan(e: MouseEvent, projectionTF: mat3): true {
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

  wheelZoom(e: WheelEvent, projectionTF: mat3): true {
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

  handleEvent(e: MouseEvent, projectionTF: mat3): boolean {
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
        viewChanged = this.wheelZoom(e as WheelEvent, projectionTF);
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

function attachCamera(canvas: HTMLCanvasElement): Camera {
  return new Camera(canvas);
}

export default attachCamera;
