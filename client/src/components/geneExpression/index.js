// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import * as d3 from "d3";
import { connect } from "react-redux";
import HistogramBrush from "../brushableHistogram";
import CirclePlus from "react-icons/lib/fa/plus-circle";
import * as globals from "../../globals";
// import ReactAutocomplete from "react-autocomplete"; /* http://emilebres.github.io/react-virtualized-checkbox/ */
import actions from "../../actions/";

@connect(state => {
  const metadata = _.get(state.controls.world, "obsAnnotations", null);
  const ranges = _.get(state.controls.world, "summary.obs", null);
  const initializeRanges = _.get(state.controls.world, "summary.obs");

  return {
    ranges,
    metadata,
    initializeRanges,
    userDefinedGenes: state.userDefinedGenes,
    world: state.controls.world,
    colorAccessor: state.controls.colorAccessor,
    allGeneNames: state.controls.allGeneNames
  };
})
class GeneExpression extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      svg: null,
      ctx: null,
      axes: null,
      dimensions: null,
      gene: null
    };
  }

  componentDidMount() {}
  handleBrushAction(selection) {
    this.props.dispatch({
      type: "continuous selection using parallel coords brushing",
      data: selection
    });
  }
  handleColorAction(key) {
    return () => {};
    const indexOfGene = 0; /* we only get one, this comes from server as needed now */

    const expressionMap = {};
    /*
      converts [{cellname: cell123, e}, ...]

      expressionMap = {
        cell123: [123, 2],
        cell789: [0, 8]
      }
    */
    // _.each(data.data.cells, cell => {
    //   /* this action is coming directly from the server */
    //   expressionMap[cell.cellname] = cell.e;
    // });
    //
    // const minExpressionCell = _.minBy(data.data.cells, cell => {
    //   return cell.e[indexOfGene];
    // });
    //
    // const maxExpressionCell = _.maxBy(data.data.cells, cell => {
    //   return cell.e[indexOfGene];
    // });
    // dispatch({
    // type: "color by expression"
    //   gene: gene,
    //   colorAccessor: gene,
    //   data,
    //   expressionMap,
    //   minExpressionCell,
    //   maxExpressionCell,
    //   indexOfGene
    // });
  }

  render() {
    return (
      <div>
        <input
          onChange={e => {
            this.setState({ gene: e.target.value });
          }}
          type="text"
        />
        <button
          style={{
            border: "none",
            background: "none",
            cursor: "pointer",
            outline: "none",
            marginTop: 20,
            padding: 0
          }}
          onClick={() => {
            this.props.dispatch(
              actions.requestGeneExpressionCountsPOST([this.state.gene])
            );
          }}
        >
          <CirclePlus
            style={{
              display: "inline-block",
              marginLeft: 7,
              color: globals.brightBlue,
              fontSize: 22
            }}
          />{" "}
          add gene{" "}
        </button>
        {this.props.world
          ? _.map(this.props.world.varDataCache, (value, key) => {
              if (key[0] === "_" /* private */) {
                return;
              } else {
                return (
                  <HistogramBrush
                    key={key}
                    field={key}
                    ranges={d3.extent(value)}
                    handleColorAction={this.handleColorAction(key).bind(this)}
                  />
                );
              }
            })
          : null}
      </div>
    );
  }
}

export default GeneExpression;

// <ReactAutocomplete
//   items={this.props.allGeneNames}
//   shouldItemRender={(item, value) =>
//     item.toLowerCase().indexOf(value.toLowerCase()) > -1
//   }
//   getItemValue={item => item}
//   renderItem={(item, highlighted) => (
//     <div
//       key={item}
//       style={{
//         backgroundColor: highlighted ? "#eee" : "transparent"
//       }}
//     >
//       {item}
//     </div>
//   )}
//   value={this.state.value}
//   onChange={e => this.setState({ value: e.target.value })}
//   onSelect={value => {
//     this.setState({ value });
//     this.props.dispatch(
//       actions.requestSingleGeneExpressionCountsPOST(
//         value
//       )
//     );
//   }}
// />

// <div>
//   {_.map(this.props.ranges, (value, key) => {
//     const isColorField = key.includes("color") || key.includes("Color");
//     if (value.range && key !== "CellName" && !isColorField) {
//       return (
//         <HistogramBrush
//           key={key}
//           metadataField={key}
//           ranges={value.range}
//         />
//       );
//     }
//   })}
//   </div>
