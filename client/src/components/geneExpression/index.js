// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import * as d3 from "d3";
import { connect } from "react-redux";
import HistogramBrush from "../brushableHistogram";
import * as globals from "../../globals";
import actions from "../../actions";
import { ErrorToastTopCenter } from "../framework/toasters";
import ExpressionButtons from "./expressionButtons";

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
      ErrorToastTopCenter.show({
        message: "That gene already exists"
      });
    } else if (userDefinedGenes.length > 15) {
      ErrorToastTopCenter.show({
        message:
          "That's too many genes, you can have at most 15 user defined genes"
      });
    } else if (!_.find(world.varAnnotations, { name: gene })) {
      ErrorToastTopCenter.show({
        message: "That doesn't appear to be a valid gene name."
      });
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
        <div
          style={{
            backgroundColor: globals.lightestGrey,
            padding: 10,
            marginTop: 30
          }}
        >
          <p style={globals.leftSidebarSectionHeading}>Custom genes</p>
          <div style={{ marginTop: 11 }} className="bp3-control-group">
            <div className="bp3-input-group">
              <span className="bp3-icon bp3-icon-symbol-cross" />
              <input
                onKeyDown={this.keyPress.bind(this)}
                onChange={e => {
                  this.setState({ gene: e.target.value });
                }}
                value={gene}
                type="text"
                className="bp3-input"
                placeholder="Add a custom gene"
                style={{ paddingRight: 94 }}
              />
            </div>
            <button
              type="button"
              className="bp3-button bp3-intent-primary"
              onClick={this.handleClick.bind(this)}
            >
              Add
            </button>
          </div>
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
        </div>
        <p
          style={Object.assign({}, globals.leftSidebarSectionHeading, {
            marginTop: 40
          })}
        >
          Differentially Expressed Genes
        </p>
        <ExpressionButtons />
        {differential.diffExp
          ? _.map(differential.diffExp, value => {
              const annotations = world.varAnnotations[value[0]];
              const { name } = annotations;
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
                  avgDiff={value[1]}
                  set1AvgExp={value[4]}
                  set2AvgExp={value[5]}
                />
              );
            })
          : null}
      </div>
    );
  }
}

export default GeneExpression;
