// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";
import _ from "lodash";
import Occupancy from "./occupancy";
import { countCategoryValues2D } from "../../util/stateManager/worldUtil";
import * as globals from "../../globals";

@connect(state => ({
  categoricalSelectionState: state.categoricalSelectionState,
  colorScale: state.colors.scale,
  colorAccessor: state.colors.colorAccessor,
  schema: _.get(state.world, "schema", null),
  world: state.world
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
    const isColorBy = metadataField === colorAccessor;
    let categories = null;
    let occupancy = null;

    if (isColorBy && schema) {
      categories = _.filter(schema.annotations.obs, {
        name: colorAccessor
      })[0].categories;
    }

    if (
      colorAccessor &&
      !isColorBy &&
      categoricalSelectionState[colorAccessor]
    ) {
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
            width: globals.leftSidebarWidth - 130,
            display: "flex",
            justifyContent: "space-between"
          }}
        >
          <label className="bp3-control bp3-checkbox">
            <input
              onChange={
                selected ? this.toggleOff.bind(this) : this.toggleOn.bind(this)
              }
              data-testclass="categorical-value-checkbox"
              data-testid={`categorical-value-checkbox-${metadataField}`}
              checked={selected}
              type="checkbox"
            />
            <span className="bp3-control-indicator" />
            {displayString}
          </label>
          <span style={{ flexShrink: 0 }}>
            {colorAccessor &&
            !isColorBy &&
            categoricalSelectionState[colorAccessor] ? (
              <Occupancy
                occupancy={occupancy.get(
                  category.categoryValues[categoryIndex]
                )}
                {...this.props}
              />
            ) : null}
          </span>
        </div>
        <span>
          <span
            data-testid={`categorical-value-count-${metadataField}-${count}`}
          >
            {count}
          </span>
          <svg
            style={{
              marginLeft: 5,
              width: 11,
              height: 11,
              backgroundColor:
                isColorBy && categories
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
