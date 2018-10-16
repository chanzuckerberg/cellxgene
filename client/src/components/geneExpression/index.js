// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import * as d3 from "d3";
import { connect } from "react-redux";
import { FaPlusCircle } from "react-icons/fa";
import HistogramBrush from "../brushableHistogram";
import * as globals from "../../globals";
// import ReactAutocomplete from "react-autocomplete"; /* http://emilebres.github.io/react-virtualized-checkbox/ */
import actions from "../../actions";

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
    colorAccessor: state.controls.colorAccessor,
    allGeneNames: state.controls.allGeneNames,
    differential: state.differential
  };
})
class GeneExpression extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      gene: ""
    };
  }

  keyPress(e) {
    if (e.keyCode === 13) {
      this.handleClick();
    }
  }

  handleClick() {
    const { world, dispatch, userDefinedGenes } = this.props;
    const { gene } = this.state;

    if (userDefinedGenes.indexOf(gene) !== -1) {
      console.log("That gene already exists");
    } else if (userDefinedGenes.length > 15) {
      console.log(
        "That's too many genes, you can have at most 15 user defined genes"
      );
    } else if (!_.find(world.varAnnotations, { name: gene })) {
      console.log("That doesn't appear to be a valid gene name.");
    } else {
      dispatch(actions.requestUserDefinedGene(gene));
      dispatch({
        type: "user defined gene",
        data: gene
      });
      this.setState({ gene: "" });
    }
  }

  render() {
    const { world, userDefinedGenes, differential } = this.props;
    const { gene } = this.state;

    return (
      <div>
        <input
          onKeyDown={this.keyPress.bind(this)}
          onChange={e => {
            this.setState({ gene: e.target.value });
          }}
          type="text"
          value={gene}
        />
        <button
          type="button"
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
          <FaPlusCircle
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
          ? _.map(userDefinedGenes, geneName => {
              const values = world.varDataCache[geneName];
              if (!values) {
                return null;
              }
              return (
                <HistogramBrush
                  key={geneName}
                  field={geneName}
                  ranges={d3.extent(values)}
                  isUserDefined
                />
              );
            })
          : null}
        {differential.diffExp ? <p> Differentially Expressed Genes </p> : null}
        {differential.diffExp
          ? _.map(differential.diffExp, value => {
              const annotations = world.varAnnotations[value[0]];
              const name = { annotations };
              const values = world.varDataCache[name];
              if (!values) {
                return null;
              }
              return (
                <HistogramBrush
                  key={name}
                  field={name}
                  ranges={d3.extent(values)}
                  isDiffExp
                />
              );
            })
          : null}
      </div>
    );
  }
}

export default GeneExpression;
