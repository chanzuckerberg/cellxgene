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
          justifyContent: "space-between"
        }}
      >
        <div
          style={{
            margin: 0,
            padding: 0,
            userSelect: "none"
          }}
        >
          <label className="bp3-control bp3-checkbox">
            <input
              onChange={
                selected ? this.toggleOff.bind(this) : this.toggleOn.bind(this)
              }
              checked={selected}
              type="checkbox"
            />
            <span className="bp3-control-indicator" />
            {value}
          </label>
        </div>
        <span>
          <span>{count}</span>
          <svg
            style={{
              marginLeft: 5,
              width: 11,
              height: 11,
              backgroundColor: c ? colorScale(value) : "inherit"
            }}
          />
        </span>
      </div>
    );
  }
}

export default CategoryValue;
