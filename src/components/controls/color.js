import React from 'react';
import _ from "lodash";
import styles from "./controls.css";
// import SectionHeader from "../framework/sectionHeader";
import { connect } from "react-redux";

const Option = ({name, selected, handleClick}) => {
  return (
    <button
      onClick={handleClick(name)}
      style={{
        textDecoration: selected ? "underline" : "none",
        border: "none",
        background: "none",
        cursor: "pointer",
        outline: "none",
      }}>
      {name}
    </button>
  )
}

@connect((state) => {
  const ranges = state.cells.cells && state.cells.cells.data.ranges ? state.cells.cells.data.ranges : null;
  const initializeRanges = state.initialize.data && state.initialize.data.data.ranges ? state.initialize.data.data.ranges : null;
  return {
    ranges,
    initializeRanges,
    colorAccessor: state.controls.colorAccessor,
  }
})
class ColorOptions extends React.Component {

  constructor(props) {
    super(props);
    this.state = {

    };
  }

  handleClick (name) {
    return () => {
      this.props.dispatch({
        type: "color changed",
        colorAccessor: name,
        rangeMaxForColorAccessor: this.props.initializeRanges[name].range.max
      });
    }
  }

  render() {
    if (!this.props.ranges || !this.props.initializeRanges) { return null }
    return (
      <div>
        Color:
        {
          _.map(this.props.ranges, (value, key) => {
            if (
              (value.options ||
              value.range) &&
              key !== "CellName"
            ) {
            return <Option
              key={key}
              handleClick={this.handleClick.bind(this)}
              selected={this.props.colorAccessor === key}
              name={key}
              />
            }
          })
        }
      </div>
    )
  }
};

export default ColorOptions;

// <SectionHeader text="Color"/>
