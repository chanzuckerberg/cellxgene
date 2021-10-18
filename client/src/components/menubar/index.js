import React from "react";
import { connect } from "react-redux";
import { ButtonGroup, AnchorButton, InputGroup, Tooltip } from "@blueprintjs/core";

import * as globals from "../../globals";
import styles from "./menubar.css";
import actions from "../../actions";
import Clip from "./clip";
import AuthButtons from "./authButtons";
import Subset from "./subset";
import UndoRedoReset from "./undoRedo";
import DiffexpButtons from "./diffexpButtons";
import Reembedding from "./reembedding";
import Preprocessing from "./preprocessing";
import { getEmbSubsetView } from "../../util/stateManager/viewStackHelpers";
import { requestSankey } from "../../actions/sankey";

@connect((state) => {
  const { annoMatrix } = state;
  const crossfilter = state.obsCrossfilter;
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
    graphInteractionMode: state.controls.graphInteractionMode,
    clipPercentileMin: Math.round(100 * (annoMatrix?.clipRange?.[0] ?? 0)),
    clipPercentileMax: Math.round(100 * (annoMatrix?.clipRange?.[1] ?? 1)),
    userDefinedGenes: state.controls.userDefinedGenes,
    colorAccessor: state.colors.colorAccessor,
    scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
    scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
    libraryVersions: state.config?.library_versions,
    auth: state.config?.authentication,
    userInfo: state.userInfo,
    undoDisabled: state["@@undoable/past"].length === 0,
    redoDisabled: state["@@undoable/future"].length === 0,
    aboutLink: state.config?.links?.["about-dataset"],
    disableDiffexp: state.config?.parameters?.["disable-diffexp"] ?? false,
    diffexpMayBeSlow:
      state.config?.parameters?.["diffexp-may-be-slow"] ?? false,
    showCentroidLabels: state.centroidLabels.showLabels,
    tosURL: state.config?.parameters?.about_legal_tos,
    privacyURL: state.config?.parameters?.about_legal_privacy,
    categoricalSelection: state.categoricalSelection,
    displaySankey: state.sankeySelection.displaySankey,
    layoutChoice: state.layoutChoice,
    outputController: state.outputController,
    sankeyController: state.sankeyController

  };
})
class MenuBar extends React.PureComponent {
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
      saveName: ""
    };
  }

  handleSankey = () => {
    const { dispatch, layoutChoice } = this.props;
    if (!layoutChoice.sankey) {
      const prom = dispatch(requestSankey());
      const links = []
      const nodes = []
      prom.then((res) => {
        let n = []
        res.edges.forEach(function (item, index) {
          links.push({
            source: item[0],
            target: item[1],
            value: res.weights[index]
          })
          n.push(item[0])
          n.push(item[1])
        });   
        n = n.filter((item, i, ar) => ar.indexOf(item) === i);
  
        n.forEach(function (item){
          nodes.push({
            id: item
          })
        })
        
        const data = {links: links, nodes: nodes}
        dispatch({type: "sankey: set data",data: data})
        dispatch({type: "toggle sankey"})
      })      
    } else {
      dispatch({type: "sankey: reset"})
      dispatch({type: "toggle sankey"})
    }

  };
  handleSaveData = () => {
    const { dispatch } = this.props;
    const { saveName } = this.state;
    dispatch(actions.requestSaveAnndataToFile(saveName))
  }
  handleReload = () => {
    const { dispatch } = this.props;
    dispatch(actions.requestReloadBackend())
  }
  isClipDisabled = () => {
    /*
    return true if clip button should be disabled.
    */
    const { pendingClipPercentiles } = this.state;
    const clipPercentileMin = pendingClipPercentiles?.clipPercentileMin;
    const clipPercentileMax = pendingClipPercentiles?.clipPercentileMax;
    const {
      clipPercentileMin: currentClipMin,
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

  handleClipOnKeyPress = (e) => {
    /*
    allow only numbers, plus other critical keys which
    may be required to make a number
    */
    if (!MenuBar.isValidDigitKeyEvent(e)) {
      e.preventDefault();
    }
  };
  handleSaveChange = (e) => {
    this.setState({saveName: e.target.value})
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
    dispatch(actions.clipAction(min, max));
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

  handleSubset = () => {
    const { dispatch } = this.props;
    dispatch(actions.subsetAction());
  };

  handleSubsetReset = () => {
    const { dispatch } = this.props;
    dispatch(actions.resetSubsetAction());
  };

  render() {
    const {
      dispatch,
      disableDiffexp,
      undoDisabled,
      redoDisabled,
      selectionTool,
      clipPercentileMin,
      clipPercentileMax,
      graphInteractionMode,
      showCentroidLabels,
      categoricalSelection,
      colorAccessor,
      subsetPossible,
      subsetResetPossible,
      userInfo,
      auth,
      displaySankey,
      layoutChoice,
      outputController,
      sankeyController
    } = this.props;
    const { pendingClipPercentiles, saveName } = this.state;
    const isColoredByCategorical = !!categoricalSelection?.[colorAccessor];
    const loading = !!outputController?.pendingFetch;
    const loadingSankey = !!sankeyController?.pendingFetch;

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
        <Reembedding />
        <Preprocessing />
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
            content={
              selectionTooltip === "select"
                ? "Lasso selection"
                : selectionTooltip
            }
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
          <Tooltip
            content="Show metadata information for hovered cells"
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
          >
            <AnchorButton
              type="button"
              icon="airplane"
              active={graphInteractionMode === "lidar"}
              onClick={() => {
                dispatch({
                  type: "change graph interaction mode",
                  data: "lidar",
                });
              }}
            />       
          </Tooltip>
          <Tooltip
            content="Display sankey plot from selected categories."
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
          >            
            <AnchorButton
                type="button"
                icon="duplicate"
                active={layoutChoice.sankey}
                disabled={!displaySankey && !layoutChoice.sankey}
                loading={loadingSankey}
                onClick={() => {
                  this.handleSankey()
                }}
              />           
            </Tooltip>
    
        </ButtonGroup>

        {disableDiffexp ? null : <DiffexpButtons />}
        <div
            style={{
              paddingBottom: "9px",
              paddingTop: "8px",
              width: "400px"
            }}>
            <InputGroup
                id="save-output-name"
                placeholder="Enter an output path to save the data to an h5ad file..."
                value={saveName}
                onChange={this.handleSaveChange}
            />            
          </div>  
          <div
            style={{
              paddingBottom: "9px",
              paddingTop: "8px"
            }}>  
          <Tooltip
            content="Write current subset into the backend."
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
          >   
          <AnchorButton
              type="button"
              icon="refresh"
              loading={loading}
              onClick={() => {
                this.handleReload()
              }}
            /> 
            </Tooltip> 
            <Tooltip
              content="Save current subset to an `.h5ad` file."
              position="bottom"
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
            >                                              
            <AnchorButton
                type="button"
                icon="floppy-disk"
                loading={loading}
                disabled={saveName===""}
                onClick={() => {
                  this.handleSaveData()
                }}
              /> 
            </Tooltip>
              
          </div>  
       
      </div>
    );
  }
}

export default MenuBar;
