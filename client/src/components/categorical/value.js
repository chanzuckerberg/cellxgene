// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";
import Occupancy from "./occupancy";
import * as globals from "../../globals";

@connect(state => ({
  categoricalSelection: state.categoricalSelection,
  colorScale: state.colors.scale,
  colorAccessor: state.colors.colorAccessor,
  schema: state.world?.schema,
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

  shouldComponentUpdate = nextProps => {
    const {
      metadataField,
      categoryIndex,
      categoricalSelection,
      colorAccessor
    } = this.props;
    const {
      categoricalSelection: newCategoricalSelection,
      colorAccessor: newcolorAccessor
    } = nextProps;

    const valueSelectionChange =
      newCategoricalSelection[metadataField].categorySelected[categoryIndex] !==
      categoricalSelection[metadataField].categorySelected[categoryIndex];

    const colorByChange = colorAccessor !== newcolorAccessor;

    return valueSelectionChange || colorByChange;
  };

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
      categoricalSelection,
      metadataField,
      categoryIndex,
      colorAccessor,
      colorScale,
      i,
      schema
    } = this.props;

    if (!categoricalSelection) return null;

    const category = categoricalSelection[metadataField];
    const selected = category.categorySelected[categoryIndex];
    const count = category.categoryCounts[categoryIndex];
    const value = category.categoryValues[categoryIndex];
    const displayString = String(
      category.categoryValues[categoryIndex]
    ).valueOf();

    /* this is the color scale, so add swatches below */
    const isColorBy = metadataField === colorAccessor;
    let categories = null;

    if (isColorBy && schema) {
      categories = schema.annotations.obsByName[colorAccessor]?.categories;
    }

    return (
      <div
        key={i}
        style={{
          display: "flex",
          alignItems: "baseline",
          justifyContent: "space-between"
        }}
        data-testclass="categorical-row"
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
              data-testclass="categorical-value-select"
              data-testid={`categorical-value-select-${metadataField}-${displayString}`}
              checked={selected}
              type="checkbox"
            />
            <span className="bp3-control-indicator" />
            <span
              data-testid={`categorical-value-${metadataField}-${displayString}`}
              data-testclass="categorical-value"
              style={{ wordBreak: "break-all" }}
            >
              {displayString}
            </span>
          </label>
          {/* color by continuous distribution histogram will go here... */}
          <span style={{ flexShrink: 0 }}>
            {colorAccessor && !isColorBy ? (
              <Occupancy category={category} {...this.props} />
            ) : null}
          </span>
        </div>
        <span>
          <span
            data-testclass="categorical-value-count"
            data-testid={`categorical-value-count-${metadataField}-${displayString}`}
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
