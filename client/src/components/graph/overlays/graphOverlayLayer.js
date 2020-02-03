import React, { PureComponent, cloneElement } from "react";
import { connect } from "react-redux";

import styles from "../graph.css";

export default
@connect(state => ({
  responsive: state.responsive
}))
class GraphOverlayLayer extends PureComponent {
  /*
    This component takes its children (assumed in the data coordinate space ([0, 1] range, origin in bottom left corner))
    and transforms itself multiple times resulting in screen space ([0, screenWidth/Height] range, origin in top left corner)

    Children are assigned in the graph component and must implement onDisplayChange()
   */
  constructor(props) {
    super(props);

    this.state = {
      display: {}
    };
  }

  matrixToTransformString = m => {
    /* 
      Translates the gl-matrix mat3 to SVG matrix transform style

                            mat3                    SVG Transform Function          
        a  c  e       
        b  d  f / [a, b, 0, c, d, 0, e, f, 1] =>  matrix(a, b, c, d, e, f) / matrix(sx, 0, 0, sy, tx, ty) / matrix(m[0] m[3] m[1] m[4] m[6] m[7])
        0  0  1      
    */
    return `matrix(${m[0]} ${m[1]} ${m[3]} ${m[4]} ${m[6]} ${m[7]})`;
  };

  reverseMatrixScaleTransformString = m => {
    return `matrix(${1 / m[0]} 0 0 ${1 / m[4]} 0 0)`;
  };

  // This is passed to all children, should be called when an overlay's display state is toggled with the overlay name and its new display state in boolean form
  onDisplayChange = (overlay, displaying) => {
    this.setState(state => {
      return { ...state, display: { ...state.display, [overlay]: displaying } };
    });
  };

  render() {
    const {
      cameraTF,
      modelTF,
      projectionTF,
      responsive,
      graphPaddingRightLeft,
      graphPaddingTop,
      children
    } = this.props;

    if (!cameraTF) return null;

    const { display } = this.state;
    const displaying = Object.values(display).some(value => value); // check to see if at least one overlay is currently displayed

    const inverseTransform = `${this.reverseMatrixScaleTransformString(
      modelTF
    )} ${this.reverseMatrixScaleTransformString(
      cameraTF
    )} ${this.reverseMatrixScaleTransformString(
      projectionTF
    )} scale(1 2) scale(1 ${1 /
      -(responsive.height - graphPaddingTop)}) scale(2 1) scale(${1 /
      (responsive.width - graphPaddingRightLeft)} 1)`;

    const newChildren = React.Children.toArray(children);

    return (
      <svg
        className={styles.graphSVG}
        width={responsive.width - graphPaddingRightLeft}
        height={responsive.height}
        pointerEvents="none"
        style={{
          zIndex: 99,
          backgroundColor: displaying ? "rgba(255, 255, 255, 0.55)" : ""
        }}
      >
        <g
          id="canvas-transformation-group-x"
          transform={`scale(${responsive.width -
            graphPaddingRightLeft} 1) scale(.5 1) translate(1 0)`}
        >
          <g
            id="canvas-transformation-group-y"
            transform={`scale(1 ${-(
              responsive.height - graphPaddingTop
            )}) translate(0 -1) scale(1 .5) translate(0 1)`}
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
                  {newChildren.map(child =>
                    cloneElement(child, { inverseTransform })
                  )}
                </g>
              </g>
            </g>
          </g>
        </g>
      </svg>
    );
  }
}
