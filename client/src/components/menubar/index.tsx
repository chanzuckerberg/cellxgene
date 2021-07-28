import React from "react";
import { connect } from "react-redux";
import { ButtonGroup, AnchorButton, Tooltip } from "@blueprintjs/core";

import * as globals from "../../globals";
// @ts-expect-error ts-migrate(2307) FIXME: Cannot find module './menubar.css' or its correspo... Remove this comment to see the full error message
import styles from "./menubar.css";
import actions from "../../actions";
import Clip from "./clip";

import AuthButtons from "./authButtons";
import Subset from "./subset";
import UndoRedoReset from "./undoRedo";
import DiffexpButtons from "./diffexpButtons";
import { getEmbSubsetView } from "../../util/stateManager/viewStackHelpers";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => {
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Defa... Remove this comment to see the full error message
  const { annoMatrix } = state;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const crossfilter = (state as any).obsCrossfilter;
  const selectedCount = crossfilter.countSelected();

  const subsetPossible =
    selectedCount !== 0 && selectedCount !== crossfilter.size(); // ie, not all and not none are selected
  const embSubsetView = getEmbSubsetView(annoMatrix);
  const subsetResetPossible = !embSubsetView
    ? annoMatrix.nObs !== annoMatrix.schema.dataframe.nObs
    : annoMatrix.nObs !== embSubsetView.nObs;

  return {
    subsetPossible,
    subsetResetPossible,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    graphInteractionMode: (state as any).controls.graphInteractionMode,
    clipPercentileMin: Math.round(100 * (annoMatrix?.clipRange?.[0] ?? 0)),
    clipPercentileMax: Math.round(100 * (annoMatrix?.clipRange?.[1] ?? 1)),
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    userDefinedGenes: (state as any).controls.userDefinedGenes,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    colorAccessor: (state as any).colors.colorAccessor,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    scatterplotXXaccessor: (state as any).controls.scatterplotXXaccessor,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    scatterplotYYaccessor: (state as any).controls.scatterplotYYaccessor,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    libraryVersions: (state as any).config?.library_versions,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    auth: (state as any).config?.authentication,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    userInfo: (state as any).userInfo,
    // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
    undoDisabled: state["@@undoable/past"].length === 0,
    // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
    redoDisabled: state["@@undoable/future"].length === 0,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    aboutLink: (state as any).config?.links?.["about-dataset"],
    disableDiffexp:
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (state as any).config?.parameters?.["disable-diffexp"] ?? false,
    diffexpMayBeSlow:
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (state as any).config?.parameters?.["diffexp-may-be-slow"] ?? false,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    showCentroidLabels: (state as any).centroidLabels.showLabels,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    tosURL: (state as any).config?.parameters?.about_legal_tos,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    privacyURL: (state as any).config?.parameters?.about_legal_privacy,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    categoricalSelection: (state as any).categoricalSelection,
  };
})
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class MenuBar extends React.PureComponent<{}, State> {
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  static isValidDigitKeyEvent(e: any) {
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

  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      pendingClipPercentiles: null,
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  isClipDisabled = () => {
    /*
    return true if clip button should be disabled.
    */
    const { pendingClipPercentiles } = this.state;
    const clipPercentileMin = pendingClipPercentiles?.clipPercentileMin;
    const clipPercentileMax = pendingClipPercentiles?.clipPercentileMax;
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'clipPercentileMin' does not exist on typ... Remove this comment to see the full error message
      clipPercentileMin: currentClipMin,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'clipPercentileMax' does not exist on typ... Remove this comment to see the full error message
      clipPercentileMax: currentClipMax,
    } = this.props;

    // if you change this test, be careful with logic around
    // comparisons between undefined / NaN handling.
    const isDisabled =
      !(clipPercentileMin < clipPercentileMax) ||
      (clipPercentileMin === currentClipMin &&
        clipPercentileMax === currentClipMax);

    return isDisabled;
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleClipOnKeyPress = (e: any) => {
    /*
    allow only numbers, plus other critical keys which
    may be required to make a number
    */
    if (!MenuBar.isValidDigitKeyEvent(e)) {
      e.preventDefault();
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleClipPercentileMinValueChange = (v: any) => {
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleClipPercentileMaxValueChange = (v: any) => {
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleClipCommit = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    const { pendingClipPercentiles } = this.state;
    const { clipPercentileMin, clipPercentileMax } = pendingClipPercentiles;
    const min = clipPercentileMin / 100;
    const max = clipPercentileMax / 100;
    dispatch(actions.clipAction(min, max));
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleClipOpening = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'clipPercentileMin' does not exist on typ... Remove this comment to see the full error message
    const { clipPercentileMin, clipPercentileMax } = this.props;
    this.setState({
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax },
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleClipClosing = () => {
    this.setState({ pendingClipPercentiles: null });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleCentroidChange = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, showCentroidLabels } = this.props;

    dispatch({
      type: "show centroid labels for category",
      showLabels: !showCentroidLabels,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleSubset = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch(actions.subsetAction());
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleSubsetReset = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch(actions.resetSubsetAction());
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
      dispatch,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'disableDiffexp' does not exist on type '... Remove this comment to see the full error message
      disableDiffexp,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'undoDisabled' does not exist on type 'Re... Remove this comment to see the full error message
      undoDisabled,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'redoDisabled' does not exist on type 'Re... Remove this comment to see the full error message
      redoDisabled,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'selectionTool' does not exist on type 'R... Remove this comment to see the full error message
      selectionTool,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'clipPercentileMin' does not exist on typ... Remove this comment to see the full error message
      clipPercentileMin,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'clipPercentileMax' does not exist on typ... Remove this comment to see the full error message
      clipPercentileMax,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'graphInteractionMode' does not exist on ... Remove this comment to see the full error message
      graphInteractionMode,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'showCentroidLabels' does not exist on ty... Remove this comment to see the full error message
      showCentroidLabels,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoricalSelection' does not exist on ... Remove this comment to see the full error message
      categoricalSelection,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorAccessor' does not exist on type 'R... Remove this comment to see the full error message
      colorAccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'subsetPossible' does not exist on type '... Remove this comment to see the full error message
      subsetPossible,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'subsetResetPossible' does not exist on t... Remove this comment to see the full error message
      subsetResetPossible,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'userInfo' does not exist on type 'Readon... Remove this comment to see the full error message
      userInfo,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'auth' does not exist on type 'Readonly<{... Remove this comment to see the full error message
      auth,
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
        <AuthButtons {...{ auth, userInfo }} />
        <UndoRedoReset
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ dispatch: any; undoDisabled: any; redoDisa... Remove this comment to see the full error message
          dispatch={dispatch}
          undoDisabled={undoDisabled}
          redoDisabled={redoDisabled}
        />
        <Clip
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ pendingClipPercentiles: any; clipPercentil... Remove this comment to see the full error message
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
              // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'IconName ... Remove this comment to see the full error message
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
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ subsetPossible: any; subsetResetPossible: ... Remove this comment to see the full error message
          subsetPossible={subsetPossible}
          subsetResetPossible={subsetResetPossible}
          handleSubset={this.handleSubset}
          handleSubsetReset={this.handleSubsetReset}
        />
        {disableDiffexp ? null : <DiffexpButtons />}
      </div>
    );
  }
}

export default MenuBar;
