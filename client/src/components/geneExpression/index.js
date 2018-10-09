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
    userDefinedGenes: state.controls.userDefinedGenes,
    world: state.controls.world,
    dimensionMap: state.controls.dimensionMap,
    colorAccessor: state.controls.colorAccessor,
    allGeneNames: state.controls.allGeneNames,
    differential: state.differential
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
      gene: ""
    };
  }
  keyPress(e) {
    if (e.keyCode == 13) {
      this.handleClick();
    }
  }
  handleClick() {
    const { world, dispatch, userDefinedGenes } = this.props;

    if (userDefinedGenes.indexOf(this.state.gene) !== -1) {
      console.log("That gene already exists");
    } else if (userDefinedGenes.length > 15) {
      console.log(
        "That's too many genes, you can have at most 15 user defined genes"
      );
    } else {
      dispatch(actions.requestGeneExpressionCountsPOST([this.state.gene]));
      dispatch({
        type: "user defined gene",
        data: this.state.gene
      });
      this.setState({ gene: "" });
    }
  }
  render() {
    const { world, userDefinedGenes } = this.props;
    return (
      <div>
        <input
          onKeyDown={this.keyPress.bind(this)}
          onChange={e => {
            this.setState({ gene: e.target.value });
          }}
          type="text"
          value={this.state.gene}
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
          onClick={this.handleClick.bind(this)}
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
        {world && userDefinedGenes.length > 0 ? (
          <p>User defined genes</p>
        ) : null}
        {world && userDefinedGenes.length > 0
          ? _.map(userDefinedGenes, (geneName, index) => {
              if (!world.varDataCache[geneName]) {
                return null;
              } else {
                const values = world.varDataCache[geneName];
                return (
                  <HistogramBrush
                    key={geneName}
                    field={geneName}
                    ranges={d3.extent(values)}
                    isUserDefined
                  />
                );
              }
            })
          : null}
        {this.props.differential.diffExp ? (
          <p> Differentially Expressed Genes </p>
        ) : null}
        {this.props.differential.diffExp
          ? _.map(this.props.differential.diffExp, (value, key) => {
              const name = world.varAnnotations[value[0]].name;
              const values = world.varDataCache[name];
              if (!world.varDataCache[world.varAnnotations[value[0]].name]) {
                return null;
              } else {
                return (
                  <HistogramBrush
                    key={name}
                    field={name}
                    ranges={d3.extent(values)}
                    isDiffExp
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
