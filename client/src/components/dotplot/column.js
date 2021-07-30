import React from "react";
import { connect, shallowEqual } from "react-redux";

import * as d3 from "d3";

import memoize from "memoize-one";

import Async from "react-async";
import ErrorLoading from "./err";
import StillLoading from "./load";
import Dot from "./dot";

import { createCategorySummaryFromDfCol } from "../../util/stateManager/controlsHelpers";

import { createColorQuery } from "../../util/stateManager/colorHelpers";

@connect((state) => ({
  annoMatrix: state.annoMatrix,
  colors: state.colors,
  genesets: state.genesets.genesets,
  pointDilation: state.pointDilation,
  differential: state.differential,
  dotplot: state.dotplot,
}))
class Column extends React.Component {
  static watchAsync(props, prevProps) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  createCategorySummaryFromDfCol = memoize(createCategorySummaryFromDfCol);

  fetchAsyncProps = async (props) => {
    const {
      annoMatrix,
      colors,
      _geneSymbol,
      _geneIndex,
      metadataField,
    } = props.watchProps;

    const [categoryData, categorySummary, colorData] = await this.fetchData(
      annoMatrix,
      metadataField,
      colors,
      _geneSymbol,
      _geneIndex
    );

    return {
      categoryData,
      categorySummary,
      colorData,
    };
  };

  async fetchData(annoMatrix, metadataField, colors, _geneSymbol) {
    /*
      fetch our data and the color-by data if appropriate, and then build a summary
      of our category and a color table for the color-by annotation.
      */
    const { schema } = annoMatrix;
    const { colorMode } = colors;
    const { genesets, differential } = this.props;
    let colorDataPromise = Promise.resolve(null);

    const query = createColorQuery(
      colorMode,
      _geneSymbol,
      schema,
      genesets,
      differential.diffExp
    );
    if (query) colorDataPromise = annoMatrix.fetch(...query);

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

  render() {
    const {
      annoMatrix,
      pointDilation,
      colors,
      viewport,
      _geneSymbol,
      _geneIndex,
      rowColumnSize,
      metadataField,
    } = this.props;

    return (
      <g key={_geneSymbol}>
        <Async
          watchFn={Column.watchAsync}
          promiseFn={this.fetchAsyncProps}
          watchProps={{
            annoMatrix,
            pointDilation,
            colors,
            viewport,
            _geneSymbol,
            _geneIndex,
            metadataField,
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
              const { categoryData, categorySummary, colorData } = asyncProps;

              if (!_geneSymbol || !colorData) return null;

              /* TODO(colinmegill) #632 wire to dotplot */

              const groupBy = categoryData.col(metadataField);
              const col = colorData.icol(0);
              const range = col.summarize();

              const histogramMap = col.histogram(
                100,
                [range.min, range.max],
                groupBy
              );

              const categories =
                annoMatrix?.schema?.annotations?.obsByName[metadataField]
                  ?.categories;
              const cellCategories = groupBy.asArray();
              const geneExpressions = col.asArray();
              let mean;
              const meanGeneExpressions = {};
              for (const c of categories) {
                const arr = [];
                for (let i = 0; i < geneExpressions.length; i += 1) {
                  if (cellCategories[i] === c) {
                    arr.push(geneExpressions[i]);
                  }
                }
                mean = arr.reduce((a, b) => a + b) / arr.length;
                meanGeneExpressions[c] = mean;
              }

              const columnColorScale = d3
                .scaleLinear()
                .domain(d3.extent(Object.values(meanGeneExpressions)))
                .range([1, 0]);

              return categorySummary.categoryValues.map(
                (val, _categoryValueIndex) => {
                  return (
                    <Dot
                      key={val}
                      categoryValue={val}
                      _categoryValueIndex={_categoryValueIndex}
                      histogramMap={histogramMap}
                      _geneSymbol={_geneSymbol}
                      _geneIndex={_geneIndex}
                      colorData={colorData}
                      rowColumnSize={rowColumnSize}
                      columnColorScale={columnColorScale}
                      meanGeneExpression={meanGeneExpressions[val]}
                    />
                  );
                }
              );
            }}
          </Async.Fulfilled>
        </Async>
      </g>
    );
  }
}

export default Column;

// updateColorTable(colors, colorDf) {
//   const { annoMatrix } = this.props;
//   const { schema } = annoMatrix;

//   /* update color table state */
//   if (!colors || !colorDf) {
//     return createColorTable(
//       null, // default mode
//       null,
//       null,
//       schema,
//       null
//     );
//   }

//   const { colorAccessor, userColors, colorMode } = colors;
//   return createColorTable(
//     colorMode,
//     colorAccessor /* TODO(colinmegill) #632 dotplot wiring */,
//     colorDf,
//     schema,
//     userColors
//   );
// }

// createColorByQuery(colors) {
//   const { annoMatrix, genesets, differential } = this.props;
//   const { schema } = annoMatrix;
//   const { colorMode, colorAccessor } = colors;

//   return createColorQuery(
//     colorMode,
//     colorAccessor /* TODO(colinmegill) #632 dotplot wiring */,
//     schema,
//     genesets,
//     differential.diffExp
//   );
// }
