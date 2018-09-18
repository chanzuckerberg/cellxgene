// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";

@connect(state => ({
  categoricalAsBooleansMap: state.controls.categoricalAsBooleansMap,
  colorScale: state.controls.colorScale,
  colorAccessor: state.controls.colorAccessor
}))
class CategoryValue extends React.Component {
  toggleOff() {
    const { dispatch, metadataField, value } = this.props;
    dispatch({
      type: "categorical metadata filter deselect",
      metadataField,
      value
    });
  }

  toggleOn() {
    const { dispatch, metadataField, value } = this.props;
    dispatch({
      type: "categorical metadata filter select",
      metadataField,
      value
    });
  }

  render() {
    const {
      categoricalAsBooleansMap,
      metadataField,
      count,
      value,
      colorAccessor,
      colorScale,
      i
    } = this.props;

    if (!categoricalAsBooleansMap) return null;

    const selected = categoricalAsBooleansMap[metadataField][value];
    /* this is the color scale, so add swatches below */
    const c = metadataField === colorAccessor;

    return (
      <div
        key={i}
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
          {value}
        </p>
        <p
          style={{
            padding: "1px 10px",
            width: 80,
            textAlign: "center",
            backgroundColor: c ? colorScale(value) : "inherit",
            color: c ? "white" : "black",
            margin: 0
            // lineHeight: "1em"
          }}
        >
          {count}
        </p>
      </div>
    );
  }
}

export default CategoryValue;
