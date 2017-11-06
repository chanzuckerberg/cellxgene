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

  return {
    ranges,
    color: state.controls.color,
  }
})
class ColorOptions extends React.Component {

  constructor(props) {
    super(props);
    this.state = {

    };
  }

  handleClick (data) {
    return () => {
      this.props.dispatch({
        type: "color changed",
        data
      });
    }
  }

  render() {
    if (!this.props.ranges) { return null }
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
              selected={this.props.color === key}
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
