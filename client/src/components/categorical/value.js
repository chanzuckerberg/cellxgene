// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";

import {
  Button,
  Tooltip,
  InputGroup,
  Menu,
  MenuItem,
  Popover,
  Position,
  PopoverInteractionKind
} from "@blueprintjs/core";
import Occupancy from "./occupancy";
import { countCategoryValues2D } from "../../util/stateManager/worldUtil";
import * as globals from "../../globals";

@connect(state => ({
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  colorScale: state.colors.scale,
  colorAccessor: state.colors.colorAccessor,
  schema: state.world?.schema,
  world: state.world
}))
class CategoryValue extends React.Component {
  constructor(props) {
    super(props);
  }

  handleDeleteValue = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "delete label",
      metadataField,
      categoryIndex
    });
  };

  handleAddCurrentSelectionToThisLabel = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "label current cell selection",
      metadataField,
      categoryIndex
    });
  };

  handleEditValue = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "label edited",
      newLabel: "foo123",
      metadataField,
      categoryIndex
    });
  };

  activateEditLabelMode = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "activate edit label mode",
      metadataField,
      categoryIndex
    });
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
      world,
      isUserAnno,
      annotations
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
            width: globals.leftSidebarWidth - 240,
            display: "flex",
            justifyContent: "flex-start"
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
              {annotations.isEditingLabelName &&
              annotations.labelEditable.category === metadataField &&
              annotations.labelEditable.label === categoryIndex
                ? null
                : displayString}
            </span>
          </label>

          {isUserAnno &&
          annotations.isEditingLabelName &&
          annotations.labelEditable.category === metadataField &&
          annotations.labelEditable.label === categoryIndex ? (
            <InputGroup
              style={{ position: "relative", top: -1 }}
              small
              defaultValue={displayString}
              rightElement={
                <Button
                  minimal
                  style={{ position: "relative", top: -1 }}
                  type="submit"
                  icon="small-tick"
                  data-testclass="submitEdit"
                  data-testid="submitEdit"
                  onClick={this.handleEditValue}
                />
              }
            />
          ) : null}
        </div>
        <span style={{ flexShrink: 0 }}>
          {colorAccessor &&
          !annotations.isEditingLabelName &&
          !isColorBy &&
          categoricalSelection[colorAccessor] ? (
            <Occupancy
              occupancy={occupancy.get(category.categoryValues[categoryIndex])}
              {...this.props}
            />
          ) : null}
        </span>
        <span>
          <span
            data-testclass="categorical-value-count"
            data-testid={`categorical-value-count-${metadataField}-${displayString}`}
          >
            {count}
          </span>

          <svg
            display={isColorBy && categories ? "auto" : "none"}
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
          {isUserAnno ? (
            <Popover
              interactionKind={PopoverInteractionKind.HOVER}
              boundary="window"
              position={Position.RIGHT}
              content={
                <Menu>
                  <MenuItem
                    icon="plus"
                    data-testclass="handleAddCurrentSelectionToThisLabel"
                    data-testid={`handleAddCurrentSelectionToThisLabel-${metadataField}`}
                    onClick={this.handleAddCurrentSelectionToThisLabel}
                    text="Add current cell selection to this label"
                  />
                  <MenuItem
                    icon="edit"
                    text="Edit this label's name"
                    data-testclass="handleEditValue"
                    data-testid={`handleEditValue-${metadataField}`}
                    onClick={this.activateEditLabelMode}
                  />
                  <MenuItem
                    icon="delete"
                    data-testclass="handleDeleteValue"
                    data-testid={`handleDeleteValue-${metadataField}`}
                    onClick={this.handleDeleteValue}
                    text="Delete this value, and reassign all cells to type 'unknown'"
                  />
                </Menu>
              }
            >
              <Button
                style={{ marginLeft: 0, position: "relative", top: -1 }}
                data-testclass="seeActions"
                data-testid={`seeActions-${metadataField}`}
                onClick={this.handleDeleteValue}
                icon="more"
                small
                minimal
              />
            </Popover>
          ) : null}
        </span>
      </div>
    );
  }
}

export default CategoryValue;
