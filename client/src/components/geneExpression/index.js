// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import * as d3 from "d3";
import { connect } from "react-redux";
import { Button, Tooltip, MenuItem } from "@blueprintjs/core";
import { Suggest } from "@blueprintjs/select";
import HistogramBrush from "../brushableHistogram";
import * as globals from "../../globals";
import actions from "../../actions";
import { postUserErrorToast } from "../framework/toasters";
import ExpressionButtons from "./expressionButtons";

function highlightText(text, query) {
  let lastIndex = 0;
  //console.log(query.split(/\s+/).filter(word => word.length>0));
  const words = query
    .split(/\s+/)
    .filter(word => word.length > 0)
    .map(escapeRegExpChars);
  if (words.length === 0) {
    return [text];
  }
  //console.log(JSON.stringify(words.join("|","gi")));
  const regexp = new RegExp(words.join("|"), "gi");
  const tokens = [];
  //console.log(text);
  while (true) {
    const match = regexp.exec(text);
    //console.log("match data=",JSON.stringify(match));
    if (!match) {
      break;
    }
    const length = match[0].length;
    const before = text.slice(lastIndex, regexp.lastIndex - length);

    if (before.length > 0) {
      tokens.push(before);
    }

    lastIndex = regexp.lastIndex;
    tokens.push(<strong key={lastIndex}>{match[0]}</strong>);
    //console.log(tokens);
  }
  const rest = text.slice(lastIndex);
  if (rest.length > 0) {
    tokens.push(rest);
  }
  return tokens;
}

function escapeRegExpChars(text) {
  return text.replace(/([.*+?^=!:${}()|\[\]\/\\])/g, "\\$1");
}

const renderGene = (gene, { handleClick, modifiers, query }) => {
  if (!modifiers.matchesPredicate) {
    return null;
  }
  const text = gene;
  return (
    <MenuItem
      active={modifiers.active}
      disabled={modifiers.disabled}
      label={"Arbitrary label!"}
      key={gene}
      onClick={g => {
        console.log("item clicked: ", g);
      }}
      text={highlightText(text, query)}
    />
  );
};

const filterGene = (query, gene) => {
  return `${gene.toLowerCase()}`.indexOf(query.toLowerCase()) >= 0;
};

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
      postUserErrorToast("That gene already exists");
    } else if (userDefinedGenes.length > 15) {
      postUserErrorToast(
        "That's too many genes, you can have at most 15 user defined genes"
      );
    } else if (!_.find(world.varAnnotations, { name: gene })) {
      postUserErrorToast("That doesn't appear to be a valid gene name.");
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
            marginTop: 30
          }}
        >
          <p
            style={Object.assign({}, globals.leftSidebarSectionHeading, {
              paddingLeft: globals.leftSidebarSectionPadding,
              margin: 0
            })}
          >
            Selected Genes
          </p>
          <div
            style={{ padding: globals.leftSidebarSectionPadding }}
            className="bp3-control-group"
          >
            <div className="bp3-input-group bp3-fill">
              <input
                onKeyDown={this.keyPress.bind(this)}
                onChange={e => {
                  this.setState({ gene: e.target.value });
                }}
                value={gene}
                type="text"
                className="bp3-input"
                placeholder="Enter a gene name"
                style={{ paddingRight: 94 }}
              />
            </div>
            <Tooltip
              content="Add a gene to see its expression levels"
              position="bottom"
            >
              <Button intent="primary" onClick={this.handleClick.bind(this)}>
                Add
              </Button>
            </Tooltip>
          </div>
          <div>
            <p> Typeahead </p>
            <Suggest
              closeOnSelect
              openOnKeyDown
              noResults={<MenuItem disabled text="No matching genes." />}
              onItemSelect={g => {
                console.log("user typed: ", g);
              }}
              inputValueRenderer={g => {
                return g;
              }}
              itemPredicate={filterGene}
              itemRenderer={renderGene}
              items={["APOD", "CD7"]}
              popoverProps={{ minimal: true }}
            />
          </div>
          {world && userDefinedGenes.length > 0
            ? _.map(userDefinedGenes, (geneName, index) => {
                const values = world.varDataCache[geneName];
                if (!values) {
                  return null;
                }
                return (
                  <HistogramBrush
                    key={geneName}
                    field={geneName}
                    zebra={index % 2 === 0}
                    ranges={d3.extent(values)}
                    isUserDefined
                  />
                );
              })
            : null}
        </div>
        <div>
          <p
            style={Object.assign({}, globals.leftSidebarSectionHeading, {
              marginTop: 40,
              paddingLeft: globals.leftSidebarSectionPadding
            })}
          >
            Differentially Expressed Genes
          </p>
          <ExpressionButtons />
          {differential.diffExp
            ? _.map(differential.diffExp, (value, index) => {
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
                    zebra={index % 2 === 0}
                    ranges={d3.extent(values)}
                    isDiffExp
                    logFoldChange={value[1]}
                    pval={value[2]}
                    pvalAdj={value[3]}
                  />
                );
              })
            : null}
        </div>
      </div>
    );
  }
}

export default GeneExpression;
