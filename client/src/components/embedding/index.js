import React from "react";
import { connect } from "react-redux";
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

@connect((state) => ({
  layoutChoice: state.layoutChoice,
  schema: state.annoMatrix?.schema,
  crossfilter: state.obsCrossfilter,
}))
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
                {schema?.dataframe?.nObs} total cells{" "}
                {/* BRUCE to extend 1559 */}
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
              <RadioGroup
                label="Embedding Choice"
                onChange={this.handleLayoutChoiceChange}
                selectedValue={layoutChoice.current}
              >
                {layoutChoice.available.map((name) => (
                  <Radio label={name} value={name} key={name} />
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
