// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import { ButtonGroup, AnchorButton, Tooltip } from "@blueprintjs/core";

import * as globals from "../../globals";
import styles from "./menubar.css";
import actions from "../../actions";
import Clip from "./clip";
import Embedding from "./embedding";
import InformationMenu from "./infoMenu";
import Subset from "./subset";
import UndoRedoReset from "./undoRedo";
import DiffexpButtons from "./diffexpButtons";

@connect((state) => ({
  universe: state.universe,
  world: state.world,
  crossfilter: state.crossfilter,
  differential: state.differential,
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
  libraryVersions: state.config?.["library_versions"],
  undoDisabled: state["@@undoable/past"].length === 0,
  redoDisabled: state["@@undoable/future"].length === 0,
  aboutLink: state.config?.links?.["about-dataset"],
  disableDiffexp: state.config?.parameters?.["disable-diffexp"] ?? false,
  diffexpMayBeSlow: state.config?.parameters?.["diffexp-may-be-slow"] ?? false,
  showCentroidLabels: state.centroidLabels.showLabels,
  tosURL: state.config?.parameters?.["about_legal_tos"],
  privacyURL: state.config?.parameters?.["about_legal_privacy"],
  categoricalSelection: state.categoricalSelection,
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
      pendingClipPercentiles: null,
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

  handleClipOnKeyPress = (e) => {
    /*
    allow only numbers, plus other critical keys which
    may be required to make a number
    */
    if (!MenuBar.isValidDigitKeyEvent(e)) {
      e.preventDefault();
    }
  };

  handleClipPercentileMinValueChange = (v) => {
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
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax },
    });
  };

  handleClipPercentileMaxValueChange = (v) => {
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
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax },
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
      clipQuantiles: { min, max },
    });
  };

  handleClipOpening = () => {
    const { clipPercentileMin, clipPercentileMax } = this.props;
    this.setState({
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax },
    });
  };

  handleClipClosing = () => {
    this.setState({ pendingClipPercentiles: null });
  };

  handleCentroidChange = () => {
    const { dispatch, showCentroidLabels } = this.props;

    dispatch({
      type: "show centroid labels for category",
      showLabels: !showCentroidLabels,
    });
  };

  subsetPossible = () => {
    const { crossfilter } = this.props;
    return (
      crossfilter.countSelected() !== 0 &&
      crossfilter.countSelected() !== crossfilter.size()
    );
  };

  subsetResetPossible = () => {
    const { world, universe } = this.props;
    return world.nObs !== universe.nObs;
  };

  render() {
    const {
      dispatch,
      libraryVersions,
      disableDiffexp,
      undoDisabled,
      redoDisabled,
      selectionTool,
      clipPercentileMin,
      clipPercentileMax,
      graphInteractionMode,
      aboutLink,
      showCentroidLabels,
      privacyURL,
      tosURL,
      categoricalSelection,
      colorAccessor,
    } = this.props;
    const { pendingClipPercentiles } = this.state;

    const isColoredByCategorical = !!categoricalSelection?.[colorAccessor];

    // constants used to create selection tool button
    const [selectionTooltip, selectionButtonIcon] =
      selectionTool === "brush"
        ? ["Brush selection", "Lasso selection"]
        : ["select", "polygon-filter"];

    return (
      <div
        style={{
          position: "absolute",
          right: 8,
          top: 0,
          display: "flex",
          flexDirection: "row-reverse",
          alignItems: "flex-start",
          flexWrap: "wrap",
          justifyContent: "flex-start",
          zIndex: 3,
        }}
      >
        <InformationMenu
          libraryVersions={libraryVersions}
          aboutLink={aboutLink}
          tosURL={tosURL}
          privacyURL={privacyURL}
        />
        <UndoRedoReset
          dispatch={dispatch}
          undoDisabled={undoDisabled}
          redoDisabled={redoDisabled}
        />
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
        <Embedding />
        <Tooltip
          content="When a category is colored by, show labels on the graph"
          position="bottom"
          disabled={graphInteractionMode === "zoom"}
        >
          <AnchorButton
            className={styles.menubarButton}
            type="button"
            data-testid="centroid-label-toggle"
            icon="property"
            onClick={this.handleCentroidChange}
            active={showCentroidLabels}
            intent={showCentroidLabels ? "primary" : "none"}
            disabled={!isColoredByCategorical}
          />
        </Tooltip>
        <ButtonGroup className={styles.menubarButton}>
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
                  data: "select",
                });
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
                  data: "zoom",
                });
              }}
            />
          </Tooltip>
        </ButtonGroup>
        <Subset
          subsetPossible={this.subsetPossible()}
          subsetResetPossible={this.subsetResetPossible()}
          handleSubset={() => {
            dispatch(actions.setWorldToSelection());
            dispatch({ type: "increment graph render counter" });
          }}
          handleSubsetReset={() => {
            dispatch(actions.resetWorldToUniverse());
            dispatch({ type: "increment graph render counter" });
          }}
        />
        {disableDiffexp ? null : <DiffexpButtons />}
      </div>
    );
  }
}

export default MenuBar;
