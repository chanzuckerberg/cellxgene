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
      schema,
      world
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
    let occupancy = null;
    let occupancyContinuous = null;

    if (isColorBy && schema) {
      categories = schema.annotations.obsByName[colorAccessor]?.categories;
    }

    if (colorAccessor && !isColorBy) {
      if (categoricalSelection[colorAccessor]) {
        const groupBy = world.obsAnnotations.col(metadataField);
        occupancy = world.obsAnnotations.col(colorAccessor).histogram(groupBy);
        // console.log(
        //   "colorAccessor:",
        //   colorAccessor,
        //   "  metadataField:",
        //   metadataField,
        //   "  value:",
        //   value
        // );
        // console.log("categoricalSelec:", categoricalSelection[colorAccessor]);
        // console.log("groupBy:", groupBy);

        // console.log("histogram for value is:", occupancy.get(value));

        console.log(occupancy);
      } else if (false /* is continuous obsannotation (n_counts) */) {
        return;
      } else {
        const groupBy = world.obsAnnotations.col(metadataField);

        occupancyContinuous = world.varData
          .col(
            colorAccessor
          ) /* this is magic and col knows what kind of data is being used for histo. Could consider separate function names to make the fork explicit*/
          .histogram(
            50,
            [-0.3585, 5.648] /* hardcoded APOD range */,
            groupBy
          ); /* Because the signature changes we really need different names for histogram to differentiate signatures  */

        console.log("histogram map:", occupancyContinuous);

        // console.log(occupancyContinuous);
        // console.log(occupancyContinuous.get("F"));
      }
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
            >
              {displayString}
            </span>
          </label>
          {/* color by continuous distribution histogram will go here... */}
          <span style={{ flexShrink: 0 }}>
            {colorAccessor &&
            !isColorBy &&
            categoricalSelection[colorAccessor] ? (
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
