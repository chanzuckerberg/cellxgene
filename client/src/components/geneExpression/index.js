// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import fuzzysort from "fuzzysort";

import { connect } from "react-redux";
import {
  MenuItem,
  Button,
  FormGroup,
  InputGroup,
  ControlGroup,
  NumericInput
} from "@blueprintjs/core";
import { Suggest } from "@blueprintjs/select";
import HistogramBrush from "../brushableHistogram";
import * as globals from "../../globals";
import actions from "../../actions";
import {
  postUserErrorToast,
  keepAroundErrorToast
} from "../framework/toasters";
import ExpressionButtons from "./expressionButtons";
import finiteExtent from "../../util/finiteExtent";

const renderGene = (fuzzySortResult, { handleClick, modifiers, query }) => {
  if (!modifiers.matchesPredicate) {
    return null;
  }
  /* the fuzzysort wraps the object with other properties, like a score */
  const geneName = fuzzySortResult.target;

  return (
    <MenuItem
      active={modifiers.active}
      disabled={modifiers.disabled}
      data-testid={`suggest-menu-item-${geneName}`}
      // Use of annotations in this way is incorrect and dataset specific.
      // See https://github.com/chanzuckerberg/cellxgene/issues/483
      // label={gene.n_counts}
      key={geneName}
      onClick={g =>
        /* this fires when user clicks a menu item */
        handleClick(g)
      }
      text={geneName}
    />
  );
};

const filterGenes = (query, genes) =>
  /* fires on load, once, and then for each character typed into the input */
  fuzzysort.go(query, genes, {
    limit: 5,
    threshold: -10000 // don't return bad results
  });

@connect(state => {
  return {
    obsAnnotations: _.get(state.world, "obsAnnotations", null),
    userDefinedGenes: state.controls.userDefinedGenes,
    userDefinedGenesLoading: state.controls.userDefinedGenesLoading,
    world: state.world,
    colorAccessor: state.colors.colorAccessor,
    differential: state.differential,
    continuousPercentileMin: state.world.continuousPercentileMin,
    continuousPercentileMax: state.world.continuousPercentileMax
  };
})
class GeneExpression extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      bulkAdd: "",
      tab: "autosuggest"
    };
  }

  placeholderGeneNames() {
    /*
    return a string containing gene name suggestions for use as a user hint.
    Eg.,    Apod, Cd74, ...
    Will return a max of 3 genes, totalling 15 characters in length.
    Randomly selects gene names.

    NOTE: the random selection means it will re-render constantly.
    */
    const { world } = this.props;
    const { varAnnotations } = world;
    const geneNames = varAnnotations.col("name").asArray();
    if (geneNames.length > 0) {
      const placeholder = [];
      let len = geneNames.length;
      const maxGeneNameCount = 3;
      const maxStrLength = 15;
      len = len < maxGeneNameCount ? len : maxGeneNameCount;
      for (let i = 0, strLen = 0; i < len && strLen < maxStrLength; i += 1) {
        const deal = Math.floor(Math.random() * geneNames.length);
        const geneName = geneNames[deal];
        placeholder.push(geneName);
        strLen += geneName.length + 2; // '2' is the length of a comma and space
      }
      placeholder.push("...");
      return placeholder.join(", ");
    }
    // default - should never happen.
    return "Apod, Cd74, ...";
  }

  handleClick(g) {
    const { world, dispatch, userDefinedGenes } = this.props;
    const gene = g.target;
    if (userDefinedGenes.indexOf(gene) !== -1) {
      postUserErrorToast("That gene already exists");
    } else if (userDefinedGenes.length > 15) {
      postUserErrorToast(
        "That's too many genes, you can have at most 15 user defined genes"
      );
    } else if (world.varAnnotations.col("name").indexOf(gene) === undefined) {
      postUserErrorToast("That doesn't appear to be a valid gene name.");
    } else {
      dispatch({ type: "single user defined gene start" });
      dispatch(actions.requestUserDefinedGene(gene));
      dispatch({ type: "single user defined gene complete" });
    }
  }

  handleBulkAddClick() {
    const { world, dispatch, userDefinedGenes } = this.props;
    const { bulkAdd } = this.state;

    /*
      test:
      Apod,,, Cd74,,    ,,,    Foo,    Bar-2,,
    */
    if (bulkAdd !== "") {
      const genes = _.pull(_.uniq(bulkAdd.split(/[ ,]+/)), "");

      dispatch({ type: "bulk user defined gene start" });
      genes.forEach(gene => {
        if (gene.length === 0) {
          keepAroundErrorToast("Must enter a gene name.");
        } else if (userDefinedGenes.indexOf(gene) !== -1) {
          keepAroundErrorToast("That gene already exists");
        } else if (
          world.varAnnotations.col("name").indexOf(gene) === undefined
        ) {
          keepAroundErrorToast(
            `${gene} doesn't appear to be a valid gene name.`
          );
        } else {
          dispatch(actions.requestUserDefinedGene(gene));
        }
      });
      dispatch({ type: "bulk user defined gene complete" });
    }

    this.setState({ bulkAdd: "" });
  }

  render() {
    const {
      world,
      userDefinedGenes,
      userDefinedGenesLoading,
      differential
    } = this.props;

    const { tab, bulkAdd } = this.state;

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
            style={{
              padding: globals.leftSidebarSectionPadding
            }}
          >
            <Button
              active={tab === "autosuggest"}
              style={{ marginRight: 5 }}
              minimal
              small
              data-testid="tab-autosuggest"
              onClick={() => {
                this.setState({ tab: "autosuggest" });
              }}
            >
              Autosuggest
            </Button>
            <Button
              active={tab === "bulkadd"}
              minimal
              small
              data-testid="tab-bulkadd"
              onClick={() => {
                this.setState({ tab: "bulkadd" });
              }}
            >
              Bulk add genes
            </Button>
          </div>

          {tab === "autosuggest" ? (
            <ControlGroup
              style={{
                paddingLeft: globals.leftSidebarSectionPadding,
                paddingBottom: globals.leftSidebarSectionPadding
              }}
            >
              <Suggest
                closeOnSelect
                openOnKeyDown
                resetOnSelect
                itemDisabled={
                  userDefinedGenesLoading ? () => true : () => false
                }
                noResults={<MenuItem disabled text="No matching genes." />}
                onItemSelect={g => {
                  /* this happens on 'enter' */
                  this.handleClick(g);
                }}
                inputProps={{ "data-testid": "gene-search" }}
                inputValueRenderer={g => {
                  return "";
                }}
                itemListPredicate={filterGenes}
                itemRenderer={renderGene.bind(this)}
                items={
                  world && world.varAnnotations
                    ? world.varAnnotations.col("name").asArray()
                    : ["No genes"]
                }
                popoverProps={{ minimal: true }}
              />
              <Button
                className="bp3-button bp3-intent-primary"
                data-testid={"add-gene"}
                loading={userDefinedGenesLoading}
              >
                Add
              </Button>
            </ControlGroup>
          ) : null}
          {tab === "bulkadd" ? (
            <div style={{ paddingLeft: globals.leftSidebarSectionPadding }}>
              <form
                onSubmit={e => {
                  e.preventDefault();
                  this.handleBulkAddClick();
                }}
              >
                <FormGroup
                  helperText="Add a list of genes (comma delimited)"
                  labelFor="text-input-bulk-add"
                >
                  <ControlGroup>
                    <InputGroup
                      onChange={e => {
                        this.setState({ bulkAdd: e.target.value });
                      }}
                      id="text-input-bulk-add"
                      data-testid="text-input-bulk-add"
                      placeholder={this.placeholderGeneNames()}
                      value={bulkAdd}
                    />
                    <Button
                      intent="primary"
                      onClick={this.handleBulkAddClick.bind(this)}
                      loading={userDefinedGenesLoading}
                    >
                      Add
                    </Button>
                  </ControlGroup>
                </FormGroup>
              </form>
            </div>
          ) : null}
          {world && userDefinedGenes.length > 0
            ? _.map(userDefinedGenes, (geneName, index) => {
                const values = world.varData.col(geneName);
                if (!values) {
                  return null;
                }
                const summary = values.summarize();
                return (
                  <HistogramBrush
                    key={geneName}
                    field={geneName}
                    zebra={index % 2 === 0}
                    ranges={summary}
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
                const name = world.varAnnotations.at(value[0], "name");
                const values = world.varData.col(name);
                if (!values) {
                  return null;
                }
                const summary = values.summarize();
                return (
                  <HistogramBrush
                    key={name}
                    field={name}
                    zebra={index % 2 === 0}
                    ranges={summary}
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
