import React, { PureComponent, cloneElement } from "react";

import styles from "../graph.css";

export default class GraphOverlayLayer extends PureComponent {
  /*
    This component takes its children (assumed in the data coordinate space ([0, 1] range, origin in bottom left corner))
    and transforms itself multiple times resulting in screen space ([0, screenWidth/Height] range, origin in top left corner)

    Children are assigned in the graph component and must implement onDisplayChange()
   */
  constructor(props) {
    super(props);
    this.state = {
      display: {},
    };
  }

  matrixToTransformString = (m) => 
    /* 
      Translates the gl-matrix mat3 to SVG matrix transform style

                            mat3                    SVG Transform Function          
        a  c  e       
        b  d  f / [a, b, 0, c, d, 0, e, f, 1] =>  matrix(a, b, c, d, e, f) / matrix(sx, 0, 0, sy, tx, ty) / matrix(m[0] m[3] m[1] m[4] m[6] m[7])
        0  0  1      
    */
     `matrix(${m[0]} ${m[1]} ${m[3]} ${m[4]} ${m[6]} ${m[7]})`
  ;

  reverseMatrixScaleTransformString = (m) => `matrix(${1 / m[0]} 0 0 ${1 / m[4]} 0 0)`;

  // This is passed to all children, should be called when an overlay's display state is toggled along with the overlay name and its new display state in boolean form
  overlaySetShowing = (overlay, displaying) => {
    this.setState((state) => ({ ...state, display: { ...state.display, [overlay]: displaying } }));
  };

  render() {
    const {
      cameraTF,
      modelTF,
      projectionTF,
      children,
      handleCanvasEvent,
      width,
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
