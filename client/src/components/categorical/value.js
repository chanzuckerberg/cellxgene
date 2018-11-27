// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";
import _ from "lodash";

@connect(state => ({
  categoricalSelectionState: state.controls.categoricalSelectionState,
  colorScale: state.controls.colorScale,
  colorAccessor: state.controls.colorAccessor,
  schema: _.get(state.controls.world, "schema", null)
}))
class CategoryValue extends React.Component {
  toggleOff() {
    const { dispatch, metadataField, optionIndex } = this.props;
    dispatch({
      type: "categorical metadata filter deselect",
      metadataField,
      optionIndex
    });
  }

  toggleOn() {
    const { dispatch, metadataField, optionIndex } = this.props;
    dispatch({
      type: "categorical metadata filter select",
      metadataField,
      optionIndex
    });
  }

  render() {
    const {
      categoricalSelectionState,
      metadataField,
      optionIndex,
      colorAccessor,
      colorScale,
      i,
      schema
    } = this.props;

    if (!categoricalSelectionState) return null;

    const category = categoricalSelectionState[metadataField];
    const selected = category.optionSelected[optionIndex];
    const count = category.optionCount[optionIndex];
    const value = category.optionValue[optionIndex];
    const displayString = String(category.optionValue[optionIndex]).valueOf();

    /* this is the color scale, so add swatches below */
    const c = metadataField === colorAccessor;
    let categories = null;

    if (c && schema) {
      categories = _.filter(schema.annotations.obs, {
        name: colorAccessor
      })[0].categories;
    }

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
            {displayString}
          </label>
        </div>
        <span>
          <span>{count}</span>
          <svg
            style={{
              marginLeft: 5,
              width: 11,
              height: 11,
              backgroundColor:
                c && categories
                  ? colorScale(categories.indexOf(value))
                  : "inherit"
            }}
          />
        </span>
      </div>
    );
  }
}

export default CategoryValue;
