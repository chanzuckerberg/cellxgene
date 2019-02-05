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
  ControlGroup
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
  const gene = fuzzySortResult.obj;
  const text = gene.name;

  return (
    <MenuItem
      active={modifiers.active}
      disabled={modifiers.disabled}
      // Use of annotations in this way is incorrect and dataset specific.
      // See https://github.com/chanzuckerberg/cellxgene/issues/483
      // label={gene.n_counts}
      key={gene.name}
      onClick={g => {
        /* this fires when user clicks a menu item */
        handleClick(g);
      }}
      text={text}
    />
  );
};

const filterGenes = (query, genes) => {
  /* fires on load, once, and then for each character typed into the input */
  return fuzzysort.go(query, genes, {
    key: "name",
    limit: 5,
    threshold: -10000 // don't return bad results
  });
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
    userDefinedGenesLoading: state.controls.userDefinedGenesLoading,
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
      bulkAdd: "",
      tab: "autosuggest"
    };
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
    } else if (!_.find(world.varAnnotations, { name: gene })) {
      postUserErrorToast("That doesn't appear to be a valid gene name.");
    } else {
      dispatch(actions.requestUserDefinedGene(gene));
      dispatch({
        type: "user defined gene",
        data: gene
      });
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

      genes.forEach(gene => {
        if (userDefinedGenes.indexOf(gene) !== -1) {
          keepAroundErrorToast("That gene already exists");
        } else if (!_.find(world.varAnnotations, { name: gene })) {
          keepAroundErrorToast(
            `${gene} doesn't appear to be a valid gene name.`
          );
        } else {
          dispatch(actions.requestUserDefinedGene(gene));
          dispatch({
            type: "user defined gene",
            data: gene
          });
        }
      });
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
                inputValueRenderer={g => {
                  return "";
                }}
                itemListPredicate={filterGenes}
                itemRenderer={renderGene.bind(this)}
                items={
                  world && world.varAnnotations
                    ? world.varAnnotations
                    : [{ name: "No genes" }]
                }
                popoverProps={{ minimal: true }}
              />
              <Button
                className="bp3-button bp3-intent-primary"
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
                      placeholder="Apod, Cd74, ..."
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
                const values = world.varDataCache[geneName];
                if (!values) {
                  return null;
                }
                return (
                  <HistogramBrush
                    key={geneName}
                    field={geneName}
                    zebra={index % 2 === 0}
                    ranges={finiteExtent(values)}
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
                    ranges={finiteExtent(values)}
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
