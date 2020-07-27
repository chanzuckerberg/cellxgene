import React from "react";
import { connect } from "react-redux";
import Async from "react-async";
import {
  ButtonGroup,
  Popover,
  Button,
  Radio,
  RadioGroup,
  Tooltip,
  Position,
} from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import { getDiscreteCellEmbeddingRowIndex } from "../../util/stateManager/viewStackHelpers";

@connect((state) => {
  return {
    layoutChoice: state.layoutChoice, // TODO: really should clean up naming, s/layout/embedding/g
    schema: state.annoMatrix?.schema,
    crossfilter: state.obsCrossfilter,
  };
})
class Embedding extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {};
  }

  handleLayoutChoiceChange = (e) => {
    const { dispatch } = this.props;
    dispatch(actions.layoutChoiceAction(e.currentTarget.value));
  };

  render() {
    const { layoutChoice, schema, crossfilter } = this.props;
    const { annoMatrix } = crossfilter;
    return (
      <ButtonGroup
        style={{
          position: "absolute",
          display: "inherit",
          left: 8,
          bottom: 8,
          zIndex: 9999,
        }}
      >
        <Popover
          target={
            <Tooltip
              content="Select embedding for visualization"
              position="top"
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
            >
              <Button
                type="button"
                data-testid="layout-choice"
                icon="heatmap"
                // minimal
                id="embedding"
                style={{
                  cursor: "pointer",
                }}
              >
                {layoutChoice?.current}: {crossfilter.countSelected()} out of{" "}
                {crossfilter.size()} cells
              </Button>
            </Tooltip>
          }
          // minimal /* removes arrow */
          position={Position.TOP_LEFT}
          content={
            <div
              style={{
                display: "flex",
                justifyContent: "flex-start",
                alignItems: "flex-start",
                flexDirection: "column",
                padding: 10,
                width: 400,
              }}
            >
              <h1>Embedding Choice</h1>
              <p style={{ fontStyle: "italic" }}>
                There are {schema?.dataframe?.nObs} cells in the entire dataset.
              </p>
              <RadioGroup
                onChange={this.handleLayoutChoiceChange}
                selectedValue={layoutChoice.current}
              >
                {layoutChoice.available.map((name) => (
                  <LayoutChoice
                    annoMatrix={annoMatrix}
                    layoutName={name}
                    key={name}
                  />
                ))}
              </RadioGroup>
            </div>
          }
        />
      </ButtonGroup>
    );
  }
}

export default Embedding;

const loadEmbeddingCounts = async ({ annoMatrix, layoutName }) => {
  const embedding = await annoMatrix.fetch("emb", layoutName);
  const discreteCellIndex = getDiscreteCellEmbeddingRowIndex(embedding);
  return { embedding, discreteCellIndex };
};

const LayoutChoice = ({ annoMatrix, layoutName }) => {
  return (
    <Async
      promiseFn={loadEmbeddingCounts}
      annoMatrix={annoMatrix}
      layoutName={layoutName}
    >
      {({ data, error, isPending }) => {
        if (error) {
          /* log, as this is unexpected */
          console.error(error);
        }
        if (error || isPending) {
          /* still loading, or errored out - just omit counts (TODO: spinner?) */
          return <Radio label={`${layoutName}`} value={layoutName} />;
        }
        if (data) {
          const { embedding, discreteCellIndex } = data;
          const isAllCells = discreteCellIndex.size() === embedding.length;
          const sizeHint = `${discreteCellIndex.size()} ${
            isAllCells ? "(all) " : ""
          }cells`;
          return (
            <Radio label={`${layoutName}: ${sizeHint}`} value={layoutName} />
          );
        }
        return null;
      }}
    </Async>
  );
};
