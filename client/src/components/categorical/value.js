// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";
import * as globals from "../../globals";
import actions from "../../actions";

@connect(state => {
  return {
    categoricalAsBooleansMap: state.controls2.categoricalAsBooleansMap,
    colorScale: state.controls2.colorScale,
    colorAccessor: state.controls2.colorAccessor
  };
})
class CategoryValue extends React.Component {
  toggleOff() {
    this.props.dispatch({
      type: "categorical metadata filter deselect",
      metadataField: this.props.metadataField,
      value: this.props.value
    });
  }

  toggleOn() {
    this.props.dispatch({
      type: "categorical metadata filter select",
      metadataField: this.props.metadataField,
      value: this.props.value
    });
  }

  render() {
    if (!this.props.categoricalAsBooleansMap) return null;

    const selected = this.props.categoricalAsBooleansMap[
      this.props.metadataField
    ][this.props.value];
    const c =
      this.props.metadataField ===
      this.props
        .colorAccessor; /* this is the color scale, so add swatches below */

    return (
      <div
        key={this.props.i}
        style={{
          display: "flex",
          alignItems: "baseline",
          justifyContent: "space-between",
          fontWeight: selected ? 700 : 400
        }}
      >
        <p
          style={{
            paddingLeft: 15,
            width: 200,
            flexShrink: 0,
            margin: 0
            // lineHeight: "1em"
          }}
        >
          <input
            style={{ position: "relative", top: 1 }}
            onChange={
              selected ? this.toggleOff.bind(this) : this.toggleOn.bind(this)
            }
            checked={selected}
            type="checkbox"
          />
          {this.props.value}
        </p>
        <p
          style={{
            padding: "1px 10px",
            width: 80,
            textAlign: "center",
            backgroundColor: c
              ? this.props.colorScale(this.props.value)
              : "inherit",
            color: c ? "white" : "black",
            margin: 0
            // lineHeight: "1em"
          }}
        >
          {this.props.count}
        </p>
      </div>
    );
  }
}

export default CategoryValue;

// toggleOnlyThis() {
//   this.props.dispatch({
//     type: "categorical metadata filter none of these",
//     metadataField: this.props.metadataField,
//     value: this.props.value
//   })
// }
// <span
//   onClick={this.toggleOnlyThis.bind(this)}
//   style={{
//     fontFamily: globals.accentFont,
//     fontSize: 10,
//     fontWeight: 100,
//     fontStyle: "italic",
//     cursor: "pointer",
//   }}>
// {"only"}
// </span>

// onClick={selected ? this.toggleOff.bind(this) : this.toggleOn.bind(this)}
// toggleOn() {
//   this.props.dispatch(
//     actions.attemptCategoricalMetadataSelection(
//       this.props.metadataField,
//       this.props.value
//     ))
// }
// toggleOff() {
//   this.props.dispatch(
//     actions.attemptCategoricalMetadataDeselection(
//       this.props.metadataField,
//       this.props.value
//     ))
// }
