import React from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  ButtonGroup,
  Tooltip,
  Slider,
  MenuItem,
} from "@blueprintjs/core";
import { Select } from "@blueprintjs/select";
import * as globals from "../../globals";
import actions from "../../actions";
import styles from "./menubar.css";

@connect((state) => ({
  reembedController: state.reembedController,
  annoMatrix: state.annoMatrix,
  numPcs: state.numPcs,
}))
class Reembedding extends React.PureComponent {
  render() {
    const { dispatch, annoMatrix, reembedController, numPcs } = this.props;
    const loading = !!reembedController?.pendingFetch;
    const disabled = annoMatrix.nObs === annoMatrix.schema.dataframe.nObs;
    const tipContent = disabled
      ? "Subset cells first, then click to recompute UMAP embedding."
      : "Click to recompute UMAP embedding on the current cell subset.";

    return (
      <ButtonGroup className={styles.menubarButton}>
        <Select
          items={
            annoMatrix.schema.layers.concat(["npcs-slider"]) ||
            [] /* this is a placeholder, could be  a subcomponent to avoid this */
          }
          filterable={false}
          itemRenderer={(d, { handleClick }) => {
            return (
              <div key={d} style={{ margin: "0 auto" }}>
                {d !== "npcs-slider" ? (
                  <MenuItem
                    data-testclass="duplicate-category-wn-option"
                    onClick={handleClick}
                    key={d}
                    text={d}
                  />
                ) : (
                  <div style={{ width: "90%", margin: "0 auto" }}>
                    <Tooltip
                      content="Select the number of PCs for reembedding."
                      position="bottom"
                      hoverOpenDelay={globals.tooltipHoverOpenDelay}
                    >
                      <Slider
                        min={5}
                        max={150}
                        stepSize={1}
                        labelStepSize={145}
                        value={numPcs.npcs}
                        onChange={(val) => {
                          dispatch({
                            type: "reembed: number of pcs update",
                            num: val,
                          });
                        }}
                      />
                    </Tooltip>
                  </div>
                )}
              </div>
            );
          }}
          onItemSelect={(d) => {
            dispatch(actions.requestReembed(d, numPcs.npcs));
          }}
        >
          <Tooltip
            content={tipContent}
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
          >
            <AnchorButton
              icon="new-object"
              disabled={disabled}
              loading={loading}
              data-testid="reembedding-options"
            />
          </Tooltip>
        </Select>
      </ButtonGroup>
    );
  }
}

export default Reembedding;
