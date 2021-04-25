import React from "react";
import { connect, shallowEqual } from "react-redux";
// import _regl from "regl";
import Async from "react-async";
import memoize from "memoize-one";

import * as globals from "../../globals";
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

  createDot = (
    metadataField,
    categoryData,
    // colorAccessor,
    colorData,
    categoryValue
    // width,
    // height
  ) => {
    /*
      Knowing that colorScale is based off continuous data,
      createHistogramBins fetches the continuous data in relation to the cells relevant to the category value.
      It then separates that data into 50 bins for drawing the mini-histogram
    */
    // const groupBy = categoryData.col(metadataField);
    const col = colorData.icol(0);
    const mean = this.average(col.asArray());

    console.log(
      categoryValue,
      mean
    ); /* 

    // const histogramMap = col.histogram(
    //   50,
    //   [range.min, range.max],
    //   groupBy
    // ); /* Because the signature changes we really need different names for histogram to differentiate signatures  */

    // const bins = histogramMap.has(categoryValue)
    //   ? histogramMap.get(categoryValue)
    //   : new Array(50).fill(0);

    // const xScale = d3.scaleLinear().domain([0, bins.length]).range([0, width]);

    // const largestBin = Math.max(...bins);

    // const yScale = d3.scaleLinear().domain([0, largestBin]).range([0, height]);

    // return {
    //   xScale,
    //   yScale,
    //   bins,
    // };
    return null;
  };

  createDots = (
    metadataField,
    categoryData,
    categorySummary,
    colorAccessor = "tissue",
    colorData,
    width,
    height
  ) => {
    if (!colorData) return null;

    return categorySummary.categoryValues.map((val) => {
      return this.createDot(
        metadataField,
        categoryData,
        colorAccessor,
        colorData,
        val,
        width,
        height
      );
    });
  };

  average = (array) => array.reduce((a, b) => a + b) / array.length;

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
              return (
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
                  <circle r="5" cx="20" cy="20" />
                  {this.createDots(
                    "tissue",
                    categoryData,
                    categorySummary,
                    colorAccessor,
                    colorData,
                    width,
                    height
                  )}
                </svg>
              );
            }}
          </Async.Fulfilled>
        </Async>
      </div>
    );
  }
}

const ErrorLoading = ({ error, width, height }) => {
  console.log(error); // log to console as this is an unepected error
  return (
    <div
      style={{
        position: "fixed",
        fontWeight: 500,
        top: height / 2,
        left: globals.leftSidebarWidth + width / 2 - 50,
      }}
    >
      <span>Failure loading dotplot</span>
    </div>
  );
};

const StillLoading = ({ width, height }) => {
  /*
    Render a busy/loading indicator
    */
  return (
    <div
      style={{
        position: "fixed",
        fontWeight: 500,
        top: height / 2,
        width,
      }}
    >
      <div
        style={{
          display: "flex",
          justifyContent: "center",
          justifyItems: "center",
          alignItems: "center",
        }}
      >
        <span style={{ fontStyle: "italic" }}>Loading dotplot</span>
      </div>
    </div>
  );
};

export default Dotplot;
