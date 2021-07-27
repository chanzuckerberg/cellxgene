import React from "react";
import { connect } from "react-redux";
import { useAsync } from "react-async";
import {
  Button,
  ButtonGroup,
  H4,
  Popover,
  Position,
  Radio,
  RadioGroup,
  Tooltip,
} from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import { getDiscreteCellEmbeddingRowIndex } from "../../util/stateManager/viewStackHelpers";

type EmbeddingState = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => {
  return {
    layoutChoice: (state as any).layoutChoice,
    schema: (state as any).annoMatrix?.schema,
    crossfilter: (state as any).obsCrossfilter,
  };
})
class Embedding extends React.PureComponent<{}, EmbeddingState> {
  constructor(props: {}) {
    super(props);
    this.state = {};
  }

  handleLayoutChoiceChange = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch(actions.layoutChoiceAction(e.currentTarget.value));
  };

  render() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'layoutChoice' does not exist on type 'Re... Remove this comment to see the full error message
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
              <H4>Embedding Choice</H4>
              <p style={{ fontStyle: "italic" }}>
                There are {schema?.dataframe?.nObs} cells in the entire dataset.
              </p>
              <EmbeddingChoices
                onChange={this.handleLayoutChoiceChange}
                annoMatrix={annoMatrix}
                layoutChoice={layoutChoice}
              />
            </div>
          }
        />
      </ButtonGroup>
    );
  }
}

export default Embedding;

const loadAllEmbeddingCounts = async ({ annoMatrix, available }: any) => {
  const embeddings = await Promise.all(
    available.map((name: any) => annoMatrix.base().fetch("emb", name))
  );
  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'name' implicitly has an 'any' type.
  return available.map((name, idx) => ({
    embeddingName: name,
    embedding: embeddings[idx],
    discreteCellIndex: getDiscreteCellEmbeddingRowIndex(embeddings[idx]),
  }));
};

const EmbeddingChoices = ({ onChange, annoMatrix, layoutChoice }: any) => {
  const { available } = layoutChoice;
  const { data, error, isPending } = useAsync({
    promiseFn: loadAllEmbeddingCounts,
    annoMatrix,
    available,
  });

  if (error) {
    /* log, as this is unexpected */
    console.error(error);
  }
  if (error || isPending) {
    /* still loading, or errored out - just omit counts (TODO: spinner?) */
    return (
      <RadioGroup onChange={onChange} selectedValue={layoutChoice.current}>
        {layoutChoice.available.map((name: any) => (
          <Radio label={`${name}`} value={name} key={name} />
        ))}
      </RadioGroup>
    );
  }
  if (data) {
    return (
      <RadioGroup onChange={onChange} selectedValue={layoutChoice.current}>
        {data.map((summary: any) => {
          const { discreteCellIndex, embeddingName } = summary;
          const sizeHint = `${discreteCellIndex.size()} cells`;
          return (
            <Radio
              label={`${embeddingName}: ${sizeHint}`}
              value={embeddingName}
              key={embeddingName}
            />
          );
        })}
      </RadioGroup>
    );
  }
  return null;
};
