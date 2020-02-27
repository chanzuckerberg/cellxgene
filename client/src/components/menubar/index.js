// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import {
  Button,
  ButtonGroup,
  AnchorButton,
  Tooltip,
  Popover,
  Position,
  RadioGroup,
  Radio,
  Icon
} from "@blueprintjs/core";
import { World } from "../../util/stateManager";
import actions from "../../actions";
import CellSetButton from "./cellSetButtons";
import InformationMenu from "./infoMenu";
import UndoRedoReset from "./undoRedoReset";
import Clip from "./clip";
import * as globals from "../../globals";

@connect(state => ({
  universe: state.universe,
  world: state.world,
  crossfilter: state.crossfilter,
  differential: state.differential,
  resettingInterface: state.controls.resettingInterface,
  layoutChoice: state.layoutChoice,
  graphInteractionMode: state.controls.graphInteractionMode,
  clipPercentileMin: Math.round(100 * (state.world?.clipQuantiles?.min ?? 0)),
  clipPercentileMax: Math.round(100 * (state.world?.clipQuantiles?.max ?? 1)),
  userDefinedGenes: state.controls.userDefinedGenes,
  diffexpGenes: state.controls.diffexpGenes,
  colorAccessor: state.colors.colorAccessor,
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
  celllist1: state.differential.celllist1,
  celllist2: state.differential.celllist2,
  libraryVersions: state.config?.library_versions, // eslint-disable-line camelcase
  undoDisabled: state["@@undoable/past"].length === 0,
  redoDisabled: state["@@undoable/future"].length === 0,
  aboutLink: state.config?.links?.["about-dataset"],
  disableDiffexp: state.config?.parameters?.["disable-diffexp"] ?? false,
  diffexpMayBeSlow: state.config?.parameters?.["diffexp-may-be-slow"] ?? false,
  showCentroidLabels: state.centroidLabels.showLabels
}))
class MenuBar extends React.Component {
  static isValidDigitKeyEvent(e) {
    /*
    Return true if this event is necessary to enter a percent number input.
    Return false if not.

    Returns true for events with keys: backspace, control, alt, meta, [0-9],
    or events that don't have a key.
    */
    if (e.key === null) return true;
    if (e.ctrlKey || e.altKey || e.metaKey) return true;

    // concept borrowed from blueprint's numericInputUtils:
    // keys that print a single character when pressed have a `key` name of
    // length 1. every other key has a longer `key` name (e.g. "Backspace",
    // "ArrowUp", "Shift"). since none of those keys can print a character
    // to the field--and since they may have important native behaviors
    // beyond printing a character--we don't want to disable their effects.
    const isSingleCharKey = e.key.length === 1;
    if (!isSingleCharKey) return true;

    const key = e.key.charCodeAt(0) - 48; /* "0" */
    return key >= 0 && key <= 9;
  }

  constructor(props) {
    super(props);
    this.state = {
      pendingClipPercentiles: null
    };
  }

  isClipDisabled = () => {
    /*
    return true if clip button should be disabled.
    */
    const { pendingClipPercentiles } = this.state;
    const clipPercentileMin = pendingClipPercentiles?.clipPercentileMin;
    const clipPercentileMax = pendingClipPercentiles?.clipPercentileMax;

    const { world } = this.props;
    const currentClipMin = 100 * world?.clipQuantiles?.min;
    const currentClipMax = 100 * world?.clipQuantiles?.max;

    // if you change this test, be careful with logic around
    // comparisons between undefined / NaN handling.
    const isDisabled =
      !(clipPercentileMin < clipPercentileMax) ||
      (clipPercentileMin === currentClipMin &&
        clipPercentileMax === currentClipMax);

    return isDisabled;
  };

  isResetDisabled = () => {
    /*
    Reset should be disabled when all of the following are true:
      * nothing is selected in the crossfilter
      * world EQ universe
      * nothing is colored by
      * there are no userDefinedGenes or diffexpGenes displayed
      * scatterplot is not displayed
      * nothing in cellset1 or cellset2
      * clip percentiles are [0,100]
    */
    const {
      crossfilter,
      world,
      universe,
      userDefinedGenes,
      diffexpGenes,
      colorAccessor,
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      celllist1,
      celllist2,
      clipPercentileMin,
      clipPercentileMax
    } = this.props;

    if (!crossfilter || !world || !universe) {
      return false;
    }
    const nothingSelected = crossfilter.countSelected() === crossfilter.size();
    const nothingColoredBy = !colorAccessor;
    const noGenes = userDefinedGenes.length === 0 && diffexpGenes.length === 0;
    const scatterNotDpl = !scatterplotXXaccessor || !scatterplotYYaccessor;
    const nothingInCellsets = !celllist1 && !celllist2;

    return (
      nothingSelected &&
      World.worldEqUniverse(world, universe) &&
      nothingColoredBy &&
      noGenes &&
      scatterNotDpl &&
      nothingInCellsets &&
      clipPercentileMax === 100 &&
      clipPercentileMin === 0
    );
  };

  resetInterface = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "interface reset started"
    });
    dispatch(actions.resetInterface());
  };

  handleClipOnKeyPress = e => {
    /*
    allow only numbers, plus other critical keys which
    may be required to make a number
    */
    if (!MenuBar.isValidDigitKeyEvent(e)) {
      e.preventDefault();
    }
  };

  handleClipPercentileMinValueChange = v => {
    /*
    Ignore anything that isn't a legit number
    */
    if (!Number.isFinite(v)) return;

    const { pendingClipPercentiles } = this.state;
    const clipPercentileMax = pendingClipPercentiles?.clipPercentileMax;

    /*
    clamp to [0, currentClipPercentileMax]
    */
    if (v <= 0) v = 0;
    if (v > 100) v = 100;
    const clipPercentileMin = Math.round(v); // paranoia
    this.setState({
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax }
    });
  };

  handleClipPercentileMaxValueChange = v => {
    /*
    Ignore anything that isn't a legit number
    */
    if (!Number.isFinite(v)) return;

    const { pendingClipPercentiles } = this.state;
    const clipPercentileMin = pendingClipPercentiles?.clipPercentileMin;

    /*
    clamp to [0, 100]
    */
    if (v < 0) v = 0;
    if (v > 100) v = 100;
    const clipPercentileMax = Math.round(v); // paranoia

    this.setState({
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax }
    });
  };

  handleClipCommit = () => {
    const { dispatch } = this.props;
    const { pendingClipPercentiles } = this.state;
    const { clipPercentileMin, clipPercentileMax } = pendingClipPercentiles;
    const min = clipPercentileMin / 100;
    const max = clipPercentileMax / 100;
    dispatch({
      type: "set clip quantiles",
      clipQuantiles: { min, max }
    });
  };

  handleClipOpening = () => {
    const { clipPercentileMin, clipPercentileMax } = this.props;
    this.setState({
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax }
    });
  };

  handleClipClosing = () => {
    this.setState({ pendingClipPercentiles: null });
  };

  handleLayoutChoiceChange = e => {
    const { dispatch } = this.props;
    dispatch({
      type: "set layout choice",
      layoutChoice: e.currentTarget.value
    });
  };

  computeDiffExp = () => {
    const { dispatch, differential } = this.props;
    if (differential.celllist1 && differential.celllist2) {
      dispatch(
        actions.requestDifferentialExpression(
          differential.celllist1,
          differential.celllist2
        )
      );
    }
  };

  clearDifferentialExpression = () => {
    const { dispatch, differential } = this.props;
    dispatch({
      type: "clear differential expression",
      diffExp: differential.diffExp
    });
    dispatch({
      type: "clear scatterplot"
    });
  };

  handleCentroidChange = () => {
    const { dispatch, showCentroidLabels } = this.props;

    dispatch({
      type: "show centroid labels for category",
      showLabels: !showCentroidLabels
    });
  };

  renderDiffExp() {
    /* diffexp-related buttons may be disabled */
    const { disableDiffexp, differential, diffexpMayBeSlow } = this.props;
    if (disableDiffexp) return null;

    const haveBothCellSets =
      !!differential.celllist1 && !!differential.celllist2;

    const slowMsg = diffexpMayBeSlow
      ? " (CAUTION: large dataset - may take longer or fail)"
      : "";
    const tipMessage = `See top 10 differentially expressed genes${slowMsg}`;

    return (
      <ButtonGroup style={{ marginRight: 10 }}>
        <CellSetButton
          {...this.props} // eslint-disable-line react/jsx-props-no-spreading
          eitherCellSetOneOrTwo={1}
        />
        <CellSetButton
          {...this.props} // eslint-disable-line react/jsx-props-no-spreading
          eitherCellSetOneOrTwo={2}
        />
        {!differential.diffExp ? (
          <Tooltip
            content={tipMessage}
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelayQuick}
          >
            <AnchorButton
              disabled={!haveBothCellSets}
              intent="primary"
              data-testid="diffexp-button"
              loading={differential.loading}
              icon="left-join"
              fill
              onClick={this.computeDiffExp}
            />
          </Tooltip>
        ) : null}

        {differential.diffExp ? (
          <Tooltip
            content="Remove differentially expressed gene list and clear cell selections"
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelayQuick}
          >
            <Button
              type="button"
              fill
              intent="warning"
              onClick={this.clearDifferentialExpression}
            >
              Clear Differential Expression
            </Button>
          </Tooltip>
        ) : null}
      </ButtonGroup>
    );
  }

  render() {
    const {
      dispatch,
      crossfilter,
      resettingInterface,
      libraryVersions,
      undoDisabled,
      redoDisabled,
      selectionTool,
      clipPercentileMin,
      clipPercentileMax,
      layoutChoice,
      graphInteractionMode,
      aboutLink,
      showCentroidLabels
    } = this.props;
    const { pendingClipPercentiles } = this.state;

    // constants used to create selection tool button
    let selectionTooltip;
    let selectionButtonIcon;
    if (selectionTool === "brush") {
      selectionTooltip = "Brush selection";
      selectionButtonIcon = "select";
    } else {
      selectionTooltip = "Lasso selection";
      selectionButtonIcon = "polygon-filter";
    }

    return (
      <div
        style={{
          position: "fixed",
          right: globals.leftSidebarWidth + 8,
          top: 8,
          display: "flex"
        }}
      >
        {this.renderDiffExp()}
        <Tooltip
          content="Show only metadata and cells which are currently selected"
          position="bottom"
          hoverOpenDelay={globals.tooltipHoverOpenDelay}
        >
          <AnchorButton
            data-testid="subset-button"
            disabled={
              crossfilter &&
              (crossfilter.countSelected() === 0 ||
                crossfilter.countSelected() === crossfilter.size())
            }
            style={{
              marginRight: 10
            }}
            onClick={() => {
              dispatch(actions.regraph());
              dispatch({ type: "increment graph render counter" });
            }}
          >
            <Icon icon="double-chevron-down" />
          </AnchorButton>
        </Tooltip>
        <ButtonGroup style={{ marginRight: "10px" }}>
          <Tooltip
            content={selectionTooltip}
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
          >
            <AnchorButton
              type="button"
              data-testid="mode-lasso"
              icon={selectionButtonIcon}
              active={graphInteractionMode === "select"}
              onClick={() => {
                dispatch({
                  type: "change graph interaction mode",
                  data: "select"
                });
              }}
              style={{
                cursor: "pointer"
              }}
            />
          </Tooltip>
          <Tooltip
            content="Drag to pan, scroll to zoom"
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
          >
            <AnchorButton
              type="button"
              data-testid="mode-pan-zoom"
              icon="zoom-in"
              active={graphInteractionMode === "zoom"}
              onClick={() => {
                dispatch({
                  type: "change graph interaction mode",
                  data: "zoom"
                });
              }}
              style={{
                cursor: "pointer"
              }}
            />
          </Tooltip>
        </ButtonGroup>
        <Tooltip
          content="When a category is colored by, show labels on the graph"
          position="bottom"
          disabled={graphInteractionMode === "zoom"}
        >
          <Button
            icon="property"
            onClick={this.handleCentroidChange}
            active={showCentroidLabels}
            intent={showCentroidLabels ? "primary" : "none"}
            style={{
              marginRight: 10
            }}
          />
        </Tooltip>
        <ButtonGroup
          style={{
            marginRight: 10
          }}
        >
          <Popover
            target={
              <Tooltip
                content="Select embedding for visualization"
                position="bottom"
                hoverOpenDelay={globals.tooltipHoverOpenDelay}
              >
                <Button
                  type="button"
                  data-testid="layout-choice"
                  icon="heatmap"
                  style={{
                    cursor: "pointer"
                  }}
                />
              </Tooltip>
            }
            position={Position.BOTTOM_RIGHT}
            content={
              <div
                style={{
                  display: "flex",
                  justifyContent: "flex-start",
                  alignItems: "flex-start",
                  flexDirection: "column",
                  padding: 10
                }}
              >
                <RadioGroup
                  label="Embedding Choice"
                  onChange={this.handleLayoutChoiceChange}
                  selectedValue={layoutChoice.current}
                >
                  {layoutChoice.available.map(name => (
                    <Radio label={name} value={name} key={name} />
                  ))}
                </RadioGroup>
              </div>
            }
          />
        </ButtonGroup>
        <Clip
          pendingClipPercentiles={pendingClipPercentiles}
          clipPercentileMin={clipPercentileMin}
          clipPercentileMax={clipPercentileMax}
          handleClipOpening={this.handleClipOpening}
          handleClipClosing={this.handleClipClosing}
          handleClipCommit={this.handleClipCommit}
          isClipDisabled={this.isClipDisabled}
          handleClipOnKeyPress={this.handleClipOnKeyPress}
          handleClipPercentileMaxValueChange={
            this.handleClipPercentileMaxValueChange
          }
          handleClipPercentileMinValueChange={
            this.handleClipPercentileMinValueChange
          }
        />
        <UndoRedoReset
          dispatch={dispatch}
          isResetDisabled={this.isResetDisabled}
          resetInterface={this.resetInterface}
          resettingInterface={resettingInterface}
          undoDisabled={undoDisabled}
          redoDisabled={redoDisabled}
        />
        <InformationMenu
          libraryVersions={libraryVersions}
          aboutLink={aboutLink}
        />
      </div>
    );
  }
}

export default MenuBar;
