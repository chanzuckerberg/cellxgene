// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";
import Occupancy from "./occupancy";
import { countCategoryValues2D } from "../../util/stateManager/worldUtil";
import * as globals from "../../globals";
import {
  Button,
  Tooltip,
  Icon,
  ControlGroup,
  InputGroup
} from "@blueprintjs/core";

@connect(state => ({
  categoricalSelection: state.categoricalSelection,
  colorScale: state.colors.scale,
  colorAccessor: state.colors.colorAccessor,
  schema: state.world?.schema,
  world: state.world
}))
class CategoryValue extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isEditing: false
    };
  }

  handleDeleteValue = () => {};

  handleEditValue = () => {
    this.setState({ isEditing: true });
  };

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

    const { isEditing } = this.state;

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

    if (isColorBy && schema) {
      categories = schema.annotations.obsByName[colorAccessor]?.categories;
    }

    if (colorAccessor && !isColorBy && categoricalSelection[colorAccessor]) {
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
          {!isEditing ? (
            <label className="bp3-control bp3-checkbox">
              <input
                onChange={
                  selected
                    ? this.toggleOff.bind(this)
                    : this.toggleOn.bind(this)
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
              <Tooltip content="Edit value" position="bottom">
                <Button
                  data-testclass="handleEditValue"
                  data-testid={`handleEditValue-${metadataField}`}
                  onClick={this.handleEditValue}
                  icon="edit"
                  small
                  minimal
                />
              </Tooltip>
            </label>
          ) : (
            <ControlGroup>
              <InputGroup small />
              <Button
                small
                icon="small-tick"
                data-testclass="editValue"
                data-testid="editValue"
                onClick={this.handleEditValue}
              />
            </ControlGroup>
          )}
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
          <Tooltip content="Delete value" position="bottom">
            <Button
              data-testclass="handleDeleteValue"
              data-testid={`handleDeleteValue-${metadataField}`}
              onClick={this.handleDeleteValue}
              icon="delete"
              small
              minimal
            />
          </Tooltip>
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
