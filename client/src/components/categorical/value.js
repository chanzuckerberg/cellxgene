// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";
import _ from "lodash";
import Occupancy from "./occupancy";
import { countCategoryValues2D } from "../../util/stateManager/worldUtil";

@connect(state => ({
  categoricalSelectionState: state.controls.categoricalSelectionState,
  colorScale: state.controls.colorScale,
  colorAccessor: state.controls.colorAccessor,
  schema: _.get(state.controls.world, "schema", null),
  world: state.controls.world
}))
class CategoryValue extends React.Component {
  toggleOff() {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "categorical metadata filter deselect",
      metadataField,
      categoryIndex
    });
  }

  toggleOn() {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "categorical metadata filter select",
      metadataField,
      categoryIndex
    });
  }

  render() {
    const {
      categoricalSelectionState,
      metadataField,
      categoryIndex,
      colorAccessor,
      colorScale,
      i,
      schema,
      world
    } = this.props;

    if (!categoricalSelectionState) return null;

    const category = categoricalSelectionState[metadataField];
    const selected = category.categorySelected[categoryIndex];
    const count = category.categoryCounts[categoryIndex];
    const value = category.categoryValues[categoryIndex];
    const displayString = String(
      category.categoryValues[categoryIndex]
    ).valueOf();

    /* this is the color scale, so add swatches below */
    const c = metadataField === colorAccessor;
    let categories = null;
    let occupancy = null;

    if (c && schema) {
      categories = _.filter(schema.annotations.obs, {
        name: colorAccessor
      })[0].categories;
    }

    if (colorAccessor && !c) {
      occupancy = countCategoryValues2D(
        metadataField,
        colorAccessor,
        world.obsAnnotations
      );
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
            userSelect: "none",
            flexShrink: 2
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
        <div>
          <span>
            {colorAccessor && !c ? (
              <Occupancy
                occupancy={occupancy.get(
                  category.categoryValues[categoryIndex]
                )}
                {...this.props}
              />
            ) : null}
          </span>
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
      </div>
    );
  }
}

export default CategoryValue;
