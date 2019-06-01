// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import _ from "lodash";
import {
  Button,
  AnchorButton,
  Tooltip,
  Popover,
  Menu,
  MenuItem,
  Position,
  NumericInput,
  Divider,
  Icon,
  RadioGroup,
  Radio
} from "@blueprintjs/core";
import Logo from "../framework/logo";
import { World } from "../../util/stateManager";
import * as globals from "../../globals";
import actions from "../../actions";
import CellSetButton from "./cellSetButtons";

@connect(state => ({
  differential: state.differential,
  world: state.world,
  loading: state.controls.loading,
  datasetTitle: state.config?.displayNames?.dataset,
  crossfilter: state.crossfilter,
  resettingInterface: state.controls.resettingInterface,
  layoutChoice: state.layoutChoice,
  graphInteractionMode: state.controls.graphInteractionMode
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
      svg: null,
      tool: null,
      container: null,
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

  computeDiffExp() {
    const { dispatch, differential } = this.props;
    if (differential.celllist1 && differential.celllist2) {
      dispatch(
        actions.requestDifferentialExpression(
          differential.celllist1,
          differential.celllist2
        )
      );
    }
  }

  clearDifferentialExpression() {
    const { dispatch, differential } = this.props;
    dispatch({
      type: "clear differential expression",
      diffExp: differential.diffExp
    });
    dispatch({
      type: "clear scatterplot"
    });
  }

  render() {
    const {
      dispatch,
      differential,
      datasetTitle,
      crossfilter,
      resettingInterface,
      libraryVersions,
      undoDisabled,
      redoDisabled,
      selectionTool,
      clipPercentileMin,
      clipPercentileMax,
      layoutChoice,
      graphInteractionMode
    } = this.props;
    const { pendingClipPercentiles } = this.state;

    const haveBothCellSets =
      !!differential.celllist1 && !!differential.celllist2;

    const clipMin =
      pendingClipPercentiles?.clipPercentileMin ?? clipPercentileMin;
    const clipMax =
      pendingClipPercentiles?.clipPercentileMax ?? clipPercentileMax;
    const activeClipClass =
      clipPercentileMin > 0 || clipPercentileMax < 100
        ? " bp3-intent-warning"
        : "";

    // constants used to create selection tool button
    let selectionTooltip;
    let selectionButtonClass;
    if (selectionTool === "brush") {
      selectionTooltip = "Brush selection";
      selectionButtonClass = "bp3-icon-select";
    } else {
      selectionTooltip = "Lasso selection";
      selectionButtonClass = "bp3-icon-polygon-filter";
    }

    return (
      <div
        style={{
          paddingTop: 10,
          paddingLeft: 10,
          paddingRight: 10,
          backgroundColor: "white",
          display: "flex",
          justifyContent: "space-between"
          // boxShadow: "0px -3px 6px 2px rgba(153,153,153,0.4)",
          // borderBottom: "2px solid red"9
        }}
      >
        <div style={{ flexShrink: 0 }}>
          <Logo size={30} />
          <span
            style={{
              fontSize: 28,
              position: "relative",
              top: -6,
              fontWeight: "bold",
              marginLeft: 5,
              color: globals.logoColor,
              userSelect: "none"
            }}
          >
            cell<span
              style={{
                position: "relative",
                top: 1,
                fontWeight: 300,
                fontSize: 24
              }}
            >
              ×
            </span>gene
          </span>
          <span
            data-testid="header"
            style={{
              fontSize: 14,
              position: "relative",
              marginLeft: 7,
              top: -8
            }}
          >
            {datasetTitle}
          </span>
        </div>
        <div
          style={{
            marginRight: 10,
            marginBottom: 10,
            flexShrink: 0
          }}
        >
          <CellSetButton {...this.props} eitherCellSetOneOrTwo={1} />
          <CellSetButton {...this.props} eitherCellSetOneOrTwo={2} />
          {!differential.diffExp ? (
            <Tooltip
              content="Add two cells selections, see the top 15 differentially expressed genes between them"
              position="bottom"
            >
              <AnchorButton
                disabled={!haveBothCellSets}
                intent="primary"
                data-testid="diffexp-button"
                loading={differential.loading}
                fill
                type="button"
                onClick={this.computeDiffExp.bind(this)}
              >
                Compute Differential Expression
              </AnchorButton>
            </Tooltip>
          ) : null}
          {differential.diffExp ? (
            <Tooltip
              content="Remove differentially expressed gene list and clear cell selections"
              position="bottom"
            >
              <Button
                type="button"
                fill
                intent="warning"
                onClick={this.clearDifferentialExpression.bind(this)}
              >
                Clear Differential Expression
              </Button>
            </Tooltip>
          ) : null}
        </div>
        <div style={{ flexShrink: 0 }}>
          <div>
            <Tooltip
              content="Show only metadata and cells which are currently selected"
              position="left"
            >
              <AnchorButton
                type="button"
                data-testid="subset-button"
                disabled={
                  crossfilter &&
                  (crossfilter.countSelected() === 0 ||
                    crossfilter.countSelected() === crossfilter.size())
                }
                style={{ marginRight: 10 }}
                onClick={() => {
                  dispatch(actions.regraph());
                  dispatch({ type: "increment graph render counter" });
                }}
              >
                subset to current selection
              </AnchorButton>
            </Tooltip>
            <div className="bp3-button-group">
              <Tooltip content={selectionTooltip} position="left">
                <Button
                  type="button"
                  data-testid="mode-lasso"
                  className={`bp3-button ${selectionButtonClass}`}
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
              <Tooltip content="Pan and zoom" position="left">
                <Button
                  type="button"
                  data-testid="mode-pan-zoom"
                  className="bp3-button bp3-icon-zoom-in"
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
            </div>
            <div
              className="bp3-button-group"
              style={{
                marginLeft: 10
              }}
            >
              <Popover
                target={
                  <Button
                    type="button"
                    data-testid="layout-choice"
                    className="bp3-button bp3-icon-heatmap"
                    style={{
                      cursor: "pointer"
                    }}
                  />
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
                      label="Layout Choice"
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
            </div>

            <div
              className="bp3-button-group"
              style={{
                marginLeft: 10
              }}
            >
              <Popover
                target={
                  <Button
                    type="button"
                    data-testid="visualization-settings"
                    className={`bp3-button bp3-icon-timeline-bar-chart ${activeClipClass}`}
                    style={{
                      cursor: "pointer"
                    }}
                  />
                }
                position={Position.BOTTOM_RIGHT}
                onOpening={this.handleClipOpening}
                onClosing={this.handleClipClosing}
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
                    <div>Clip all continuous values to percentile range</div>
                    <div
                      style={{
                        display: "flex",
                        justifyContent: "space-between",
                        alignItems: "center",
                        paddingTop: 5,
                        paddingBottom: 5
                      }}
                    >
                      <NumericInput
                        style={{ width: 50 }}
                        data-testid="clip-min-input"
                        onValueChange={this.handleClipPercentileMinValueChange}
                        onKeyPress={this.handleClipOnKeyPress}
                        value={clipMin}
                        min={0}
                        max={100}
                        fill={false}
                        minorStepSize={null}
                        rightElement={
                          <div style={{ padding: "4px 2px" }}>
                            <Icon
                              icon="percentage"
                              intent="primary"
                              iconSize={14}
                            />
                          </div>
                        }
                      />
                      <span style={{ marginRight: 5, marginLeft: 5 }}> - </span>
                      <NumericInput
                        style={{ width: 50 }}
                        data-testid="clip-max-input"
                        onValueChange={this.handleClipPercentileMaxValueChange}
                        onKeyPress={this.handleClipOnKeyPress}
                        value={clipMax}
                        min={0}
                        max={100}
                        fill={false}
                        minorStepSize={null}
                        rightElement={
                          <div style={{ padding: "4px 2px" }}>
                            <Icon
                              icon="percentage"
                              intent="primary"
                              iconSize={14}
                            />
                          </div>
                        }
                      />
                      <Button
                        type="button"
                        data-testid="clip-commit"
                        className="bp3-button"
                        disabled={this.isClipDisabled()}
                        style={{
                          cursor: "pointer",
                          marginRight: 5,
                          marginLeft: 5
                        }}
                        onClick={this.handleClipCommit}
                      >
                        Clip
                      </Button>
                    </div>
                  </div>
                }
              />
            </div>
            <div
              className="bp3-button-group"
              style={{
                marginLeft: 10
              }}
            >
              <Tooltip content="Undo" position="left">
                <AnchorButton
                  type="button"
                  className="bp3-button bp3-icon-undo"
                  disabled={undoDisabled}
                  onClick={() => {
                    dispatch({ type: "@@undoable/undo" });
                  }}
                  style={{
                    cursor: "pointer"
                  }}
                />
              </Tooltip>
              <Tooltip content="Redo" position="left">
                <AnchorButton
                  type="button"
                  className="bp3-button bp3-icon-redo"
                  disabled={redoDisabled}
                  onClick={() => {
                    dispatch({ type: "@@undoable/redo" });
                  }}
                  style={{
                    cursor: "pointer"
                  }}
                />
              </Tooltip>
            </div>
            <Tooltip
              content="Reset cellxgene, clearing all selections"
              position="left"
            >
              <AnchorButton
                disabled={this.isResetDisabled()}
                style={{ marginLeft: 10 }}
                type="button"
                loading={resettingInterface}
                intent="none"
                icon="refresh"
                onClick={this.resetInterface}
                data-testid="reset"
                data-testclass={`resetting-${resettingInterface}`}
              >
                reset
              </AnchorButton>
            </Tooltip>

            <div style={{ marginLeft: 10 }} className="bp3-button-group">
              <Popover
                content={
                  <Menu>
                    <MenuItem
                      href="https://chanzuckerberg.github.io/cellxgene/faq.html"
                      target="_blank"
                      icon="help"
                      text="FAQ"
                    />
                    <MenuItem
                      href="https://join-cellxgene-users.herokuapp.com/"
                      target="_blank"
                      icon="chat"
                      text="Chat"
                    />
                    <MenuItem
                      href="https://chanzuckerberg.github.io/cellxgene/"
                      target="_blank"
                      icon="book"
                      text="Docs"
                    />
                    <MenuItem
                      href="https://github.com/chanzuckerberg/cellxgene"
                      target="_blank"
                      icon="git-branch"
                      text="Github"
                    />
                    <MenuItem
                      target="_blank"
                      text={`cellxgene v${
                        libraryVersions && libraryVersions.cellxgene
                          ? libraryVersions.cellxgene
                          : null
                      }`}
                    />
                    <MenuItem text="MIT License" />
                  </Menu>
                }
                position={Position.BOTTOM_RIGHT}
              >
                <Button
                  type="button"
                  className="bp3-button bp3-icon-info-sign"
                  style={{
                    cursor: "pointer"
                  }}
                />
              </Popover>
            </div>
          </div>
        </div>
      </div>
    );
  }
}

export default MenuBar;
