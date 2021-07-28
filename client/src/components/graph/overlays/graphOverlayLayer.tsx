import React, { PureComponent, cloneElement } from "react";

// @ts-expect-error ts-migrate(2307) FIXME: Cannot find module '../graph.css' or its correspon... Remove this comment to see the full error message
import styles from "../graph.css";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
export default class GraphOverlayLayer extends PureComponent<{}, State> {
  /*
    This component takes its children (assumed in the data coordinate space ([0, 1] range, origin in bottom left corner))
    and transforms itself multiple times resulting in screen space ([0, screenWidth/Height] range, origin in top left corner)

    Children are assigned in the graph component and must implement onDisplayChange()
   */
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      display: {},
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  matrixToTransformString = (m: any) => {
    /* 
      Translates the gl-matrix mat3 to SVG matrix transform style

                            mat3                    SVG Transform Function          
        a  c  e       
        b  d  f / [a, b, 0, c, d, 0, e, f, 1] =>  matrix(a, b, c, d, e, f) / matrix(sx, 0, 0, sy, tx, ty) / matrix(m[0] m[3] m[1] m[4] m[6] m[7])
        0  0  1      
    */
    return `matrix(${m[0]} ${m[1]} ${m[3]} ${m[4]} ${m[6]} ${m[7]})`;
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  reverseMatrixScaleTransformString = (m: any) => {
    return `matrix(${1 / m[0]} 0 0 ${1 / m[4]} 0 0)`;
  };

  // This is passed to all children, should be called when an overlay's display state is toggled along with the overlay name and its new display state in boolean form
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  overlaySetShowing = (overlay: any, displaying: any) => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    this.setState((state: any) => {
      return { ...state, display: { ...state.display, [overlay]: displaying } };
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'cameraTF' does not exist on type 'Readon... Remove this comment to see the full error message
      cameraTF,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'modelTF' does not exist on type 'Readonl... Remove this comment to see the full error message
      modelTF,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'projectionTF' does not exist on type 'Re... Remove this comment to see the full error message
      projectionTF,
      children,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'handleCanvasEvent' does not exist on typ... Remove this comment to see the full error message
      handleCanvasEvent,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'width' does not exist on type 'Readonly<... Remove this comment to see the full error message
      width,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'height' does not exist on type 'Readonly... Remove this comment to see the full error message
      height,
    } = this.props;
    const { display } = this.state;

    if (!cameraTF) return null;

    const displaying = Object.values(display).some((value) => value); // check to see if at least one overlay is currently displayed

    const inverseTransform = `${this.reverseMatrixScaleTransformString(
      modelTF
    )} ${this.reverseMatrixScaleTransformString(
      cameraTF
    )} ${this.reverseMatrixScaleTransformString(
      projectionTF
    )} scale(1 2) scale(1 ${1 / -height}) scale(2 1) scale(${1 / width} 1)`;

    // Copy the children passed with the overlay and add the inverse transform and onDisplayChange props
    const newChildren = React.Children.map(children, (child) =>
      // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
      cloneElement(child, {
        inverseTransform,
        overlaySetShowing: this.overlaySetShowing,
      })
    );

    return (
      <svg
        className={styles.graphSVG}
        width={width}
        height={height}
        pointerEvents="none"
        style={{
          position: "absolute",
          top: 0,
          left: 0,
          zIndex: 2,
          backgroundColor: displaying ? "rgba(255, 255, 255, 0.55)" : "",
        }}
        onMouseMove={handleCanvasEvent}
        onWheel={handleCanvasEvent}
      >
        <g
          id="canvas-transformation-group-x"
          transform={`scale(${width} 1) scale(.5 1) translate(1 0)`}
        >
          <g
            id="canvas-transformation-group-y"
            transform={`scale(1 ${-height}) translate(0 -1) scale(1 .5) translate(0 1)`}
          >
            <g
              id="projection-transformation-group"
              transform={this.matrixToTransformString(projectionTF)}
            >
              <g
                id="camera-transformation-group"
                transform={this.matrixToTransformString(cameraTF)}
              >
                <g
                  id="model-transformation-group"
                  transform={this.matrixToTransformString(modelTF)}
                >
                  {newChildren}
                </g>
              </g>
            </g>
          </g>
        </g>
      </svg>
    );
  }
}
