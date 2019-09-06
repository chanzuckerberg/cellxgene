// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";
import Occupancy from "./occupancy";
import * as globals from "../../globals";
import styles from "./categorical.css";
import { Tooltip } from "@blueprintjs/core";

@connect(state => ({
  categoricalSelection: state.categoricalSelection,
  colorScale: state.colors.scale,
  colorAccessor: state.colors.colorAccessor,
  schema: state.world?.schema,
  world: state.world
}))
class CategoryValue extends React.Component {
  toggleOff = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "categorical metadata filter deselect",
      metadataField,
      categoryIndex
    });
  };

  shouldComponentUpdate = nextProps => {
    /* 
    Checks to see if at least one of the following changed:
    * world state
    * the color accessor (what is currently being colored by)
    * if this catagorical value's selection status has changed
    
    If and only if true, update the component
    */
    const { props } = this;
    const { metadataField, categoryIndex, categoricalSelection } = props;
    const { categoricalSelection: newCategoricalSelection } = nextProps;

    const valueSelectionChange =
      categoricalSelection[metadataField].categoryValueSelected[
        categoryIndex
      ] !==
      newCategoricalSelection[metadataField].categoryValueSelected[
        categoryIndex
      ];

    const worldChange = props.world !== nextProps.world;
    const colorAccessorChange = props.colorAccessor !== nextProps.colorAccessor;

    return valueSelectionChange || worldChange || colorAccessorChange;
  };

  toggleOn = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "categorical metadata filter select",
      metadataField,
      categoryIndex
    });
  };

  handleMouseEnter = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "category value mouse hover start",
      metadataField,
      categoryIndex
    });
  };

  handleMouseExit = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "category value mouse hover end",
      metadataField,
      categoryIndex
    });
  };

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
    const selected = category.categoryValueSelected[categoryIndex];
    const count = category.categoryValueCounts[categoryIndex];
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

    let truncatedString = null;

    if (colorAccessor && displayString.length > 17) {
      truncatedString = `${displayString.slice(0, 8)}…${displayString.slice(
        -8
      )}`;
    } else if (displayString.length > 37) {
      truncatedString = `${displayString.slice(0, 18)}…${displayString.slice(
        -18
      )}`;
    }

    if (truncatedString) {
      console.log("TRUNCATED:", displayString);
    }

    return (
      <div
        key={i}
        className={styles.value}
        data-testclass="categorical-row"
        style={{
          padding: "4px 7px",
          display: "flex",
          alignItems: "baseline",
          justifyContent: "space-between",
          marginBottom: "2px",
          borderRadius: "2px"
        }}
        onMouseEnter={this.handleMouseEnter}
        onMouseLeave={this.handleMouseExit}
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
          <div style={{ display: "flex", alignItems: "baseline" }}>
            <label className="bp3-control bp3-checkbox" style={{ margin: 0 }}>
              <input
                onChange={selected ? this.toggleOff : this.toggleOn}
                data-testclass="categorical-value-select"
                data-testid={`categorical-value-select-${metadataField}-${displayString}`}
                checked={selected}
                type="checkbox"
              />
              <span
                className="bp3-control-indicator"
                onMouseEnter={this.handleMouseExit}
                onMouseLeave={this.handleMouseEnter}
              />
            </label>
            <Tooltip
              content={displayString}
              disabled={truncatedString === null}
            >
              <span
                data-testid={`categorical-value-${metadataField}-${displayString}`}
                data-testclass="categorical-value"
                style={{
                  overflow: "hidden",
                  height: "1.1em",
                  lineHeight: "1.1em",
                  wordBreak: "break-all",
                  display: "inline-block"
                }}
              >
                {truncatedString || displayString}
              </span>
            </Tooltip>
          </div>
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
