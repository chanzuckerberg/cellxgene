import React from "react";
import { connect, shallowEqual } from "react-redux";
// import _regl from "regl";
import Async from "react-async";
import memoize from "memoize-one";
// import * as d3 from "d3";

import ErrorLoading from "./err";
import StillLoading from "./load";

import { createCategorySummaryFromDfCol } from "../../util/stateManager/controlsHelpers";

import {
  createColorTable,
  createColorQuery,
} from "../../util/stateManager/colorHelpers";

@connect((state) => ({
  annoMatrix: state.annoMatrix,
  layoutChoice: state.layoutChoice,
  colors: state.colors,
  genesets: state.genesets.genesets,
  differential: state.differential,
  pointDilation: state.pointDilation,
}))
class Dotplot extends React.Component {
  static watchAsync(props, prevProps) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  createCategorySummaryFromDfCol = memoize(createCategorySummaryFromDfCol);

  constructor(props) {
    super(props);
    const viewport = this.getViewportDimensions();
    this.dotplotTopPadding = 150;
    this.dotplotLeftPadding = 200;

    this.state = {
      viewport,
    };
  }

  componentDidMount() {
    window.addEventListener("resize", this.handleResize);
  }

  componentWillUnmount() {
    window.removeEventListener("resize", this.handleResize);
  }

  //   handleResize = () => {};

  getViewportDimensions = () => {
    const { viewportRef } = this.props;
    return {
      height: viewportRef.clientHeight,
      width: viewportRef.clientWidth,
    };
  };

  createRow = (
    metadataField,
    categoryData,
    colorAccessor,
    colorData,
    categoryValue,
    width,
    height,
    index,
    histogramMap
  ) => {
    const bins = histogramMap.has(categoryValue)
      ? histogramMap.get(categoryValue)
      : new Array(100).fill(0);

    const zeros = bins.shift(); /* MUTATES, REMOVES FIRST ELEMENT */
    const rest = bins.reduce(
      (acc, current) => acc + current
    ); /* SUM REMAINING ELEMENTS */

    return (
      <g
        key={`${index}_${categoryValue}`}
        transform={`translate(${this.dotplotLeftPadding}, ${
          index * 15 + this.dotplotTopPadding
        })`}
      >
        <text textAnchor="end">{categoryValue}</text>
        <circle
          r={(rest / zeros) * 10}
          cx="11"
          cy="-4"
          style={{ fill: "rgba(70,130,180,.3)", stroke: "none" }}
        />
      </g>
    );
  };

  dotplot = (
    metadataField,
    categoryData,
    categorySummary,
    colorAccessor = "tissue",
    colorData,
    width,
    height
  ) => {
    if (!colorAccessor || !colorData) return null;

    const groupBy = categoryData.col(metadataField);
    const col = colorData.icol(0);
    const range = col.summarize();

    const histogramMap = col.histogram(100, [range.min, range.max], groupBy);

    return categorySummary.categoryValues.map((val, index) => {
      return this.createRow(
        metadataField,
        categoryData,
        colorAccessor,
        colorData,
        val,
        width,
        height,
        index,
        histogramMap
      );
    });
  };

  fetchAsyncProps = async (props) => {
    const { viewport, annoMatrix, colors } = props.watchProps;

    const [categoryData, categorySummary, colorData] = await this.fetchData(
      annoMatrix,
      "tissue",
      colors
    );

    const { width, height } = viewport;
    return {
      width,
      height,
      categoryData,
      categorySummary,
      colorData,
    };
  };

  async fetchData(annoMatrix, metadataField, colors) {
    /*
    fetch our data and the color-by data if appropriate, and then build a summary
    of our category and a color table for the color-by annotation.
    */
    const { schema } = annoMatrix;
    const { colorAccessor, colorMode } = colors;
    const { genesets, differential } = this.props;
    let colorDataPromise = Promise.resolve(null);
    if (colorAccessor) {
      const query = createColorQuery(
        colorMode,
        colorAccessor,
        schema,
        genesets,
        differential.diffExp
      );
      if (query) colorDataPromise = annoMatrix.fetch(...query);
    }
    const [categoryData, colorData] = await Promise.all([
      annoMatrix.fetch("obs", metadataField),
      colorDataPromise,
    ]);

    // our data
    const column = categoryData.icol(0);
    const colSchema = schema.annotations.obsByName[metadataField];

    const categorySummary = this.createCategorySummaryFromDfCol(
      column,
      colSchema
    );
    return [categoryData, categorySummary, colorData];
  }

  updateColorTable(colors, colorDf) {
    const { annoMatrix } = this.props;
    const { schema } = annoMatrix;

    /* update color table state */
    if (!colors || !colorDf) {
      return createColorTable(
        null, // default mode
        null,
        null,
        schema,
        null
      );
    }

    const { colorAccessor, userColors, colorMode } = colors;
    return createColorTable(
      colorMode,
      colorAccessor,
      colorDf,
      schema,
      userColors
    );
  }

  createColorByQuery(colors) {
    const { annoMatrix, genesets, differential } = this.props;
    const { schema } = annoMatrix;
    const { colorMode, colorAccessor } = colors;

    return createColorQuery(
      colorMode,
      colorAccessor,
      schema,
      genesets,
      differential.diffExp
    );
  }

  render() {
    const { viewport } = this.state;
    const { annoMatrix, colors, pointDilation } = this.props;

    return (
      <div
        id="dotplot-wrapper"
        style={{
          position: "relative",
          top: 0,
          left: 0,
        }}
      >
        <svg
          style={{
            position: "absolute",
            top: 0,
            left: 0,
            zIndex: 1,
          }}
          width={viewport.width}
          height={viewport.height}
        >
          <Async
            watchFn={Dotplot.watchAsync}
            promiseFn={this.fetchAsyncProps}
            watchProps={{
              annoMatrix,
              pointDilation,
              colors,
              viewport,
            }}
          >
            <Async.Pending initial>
              <StillLoading width={viewport.width} height={viewport.height} />
            </Async.Pending>
            <Async.Rejected>
              {(error) => (
                <ErrorLoading
                  width={viewport.width}
                  height={viewport.height}
                  error={error}
                />
              )}
            </Async.Rejected>
            <Async.Fulfilled persist>
              {(asyncProps) => {
                const {
                  colorAccessor,
                  categoryData,
                  width,
                  height,
                  categorySummary,
                  colorData,
                } = asyncProps;
                return this.dotplot(
                  "tissue",
                  categoryData,
                  categorySummary,
                  colorAccessor,
                  colorData,
                  width,
                  height
                );
              }}
            </Async.Fulfilled>
          </Async>
        </svg>
      </div>
    );
  }
}

export default Dotplot;
