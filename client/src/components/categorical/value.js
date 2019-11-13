// jshint esversion: 6
import { connect } from "react-redux";
import React from "react";

import {
  Button,
  InputGroup,
  Menu,
  MenuItem,
  Popover,
  Position,
  Icon,
  PopoverInteractionKind
} from "@blueprintjs/core";
import Occupancy from "./occupancy";
import * as globals from "../../globals";
import styles from "./categorical.css";
import { Tooltip } from "@blueprintjs/core";

import { AnnotationsHelpers } from "../../util/stateManager";

@connect(state => ({
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  colorScale: state.colors.scale,
  colorAccessor: state.colors.colorAccessor,
  schema: state.world?.schema,
  world: state.world,
  crossfilter: state.crossfilter
}))
class CategoryValue extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      editedLabelText: ""
    };
  }

  handleDeleteValue = () => {
    const {
      dispatch,
      metadataField,
      categoryIndex,
      categoricalSelection
    } = this.props;
    const category = categoricalSelection[metadataField];
    const label = category.categoryValues[categoryIndex];
    dispatch({
      type: "annotation: delete label",
      metadataField,
      label
    });
  };

  handleAddCurrentSelectionToThisLabel = () => {
    const {
      dispatch,
      metadataField,
      categoryIndex,
      categoricalSelection
    } = this.props;
    const category = categoricalSelection[metadataField];
    const label = category.categoryValues[categoryIndex];
    dispatch({
      type: "annotation: label current cell selection",
      metadataField,
      categoryIndex,
      label
    });
  };

  handleEditValue = () => {
    const {
      dispatch,
      metadataField,
      categoryIndex,
      categoricalSelection
    } = this.props;
    const { editedLabelText } = this.state;
    const category = categoricalSelection[metadataField];
    const label = category.categoryValues[categoryIndex];
    dispatch({
      type: "annotation: label edited",
      editedLabel: editedLabelText,
      metadataField,
      categoryIndex,
      label
    });
    this.setState({ editedLabelText: "" });
  };

  activateEditLabelMode = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "annotation: activate edit label mode",
      metadataField,
      categoryIndex
    });
  };

  cancelEdit = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "annotation: cancel edit label mode",
      metadataField,
      categoryIndex
    });
  };

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
    * the crossfilter (ie, global selection state)

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
    const annotationsChange = props.annotations !== nextProps.annotations;
    const crossfilterChange = props.crossfilter !== nextProps.crossfilter;

    return (
      valueSelectionChange ||
      worldChange ||
      colorAccessorChange ||
      annotationsChange ||
      crossfilterChange
    );
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

  isAddCurrentSelectionDisabled(category, value) {
    /*
    disable "add current selection to label", if one of the following is true:
    1. no cells are selected
    2. all currently selected cells already have this label, on this category
    */
    const { crossfilter, world } = this.props;

    // 1. no cells selected?
    if (crossfilter.countSelected() === 0) {
      return true;
    }
    // 2. all selected cells already have the label
    const mask = crossfilter.allSelectedMask();
    if (
      AnnotationsHelpers.allHaveLabelByMask(
        world.obsAnnotations,
        category,
        value,
        mask
      )
    ) {
      return true;
    }
    // else, don't disable
    return false;
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
      isUserAnno,
      annotations,
      // flippedProps is potentially brittle, their docs want {...flippedProps} on our div,
      // our lint doesn't like jsx spread, we are version pinned to prevent api change on their part
      flippedProps
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

    if (
      colorAccessor &&
      !isColorBy &&
      displayString.length > globals.categoryLabelDisplayStringShortLength
    ) {
      truncatedString = `${displayString.slice(
        0,
        globals.categoryLabelDisplayStringShortLength / 2
      )}…${displayString.slice(
        -globals.categoryLabelDisplayStringShortLength / 2
      )}`;
    } else if (
      displayString.length > globals.categoryLabelDisplayStringLongLength
    ) {
      truncatedString = `${displayString.slice(
        0,
        globals.categoryLabelDisplayStringLongLength / 2
      )}…${displayString.slice(
        -globals.categoryLabelDisplayStringLongLength / 2
      )}`;
    }

    return (
      <div
        key={i}
        data-flip-config={flippedProps["data-flip-config"]}
        data-flip-id={flippedProps["data-flip-id"]}
        data-portal-key={flippedProps["data-portal-key"]}
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
              hoverOpenDelay={globals.tooltipHoverOpenDelayQuick}
            >
              <span
                data-testid={`categorical-value-${metadataField}-${displayString}`}
                data-testclass="categorical-value"
                style={{
                  color:
                    displayString === globals.unassignedCategoryLabel
                      ? "#ababab"
                      : "black",
                  fontStyle:
                    displayString === globals.unassignedCategoryLabel
                      ? "italic"
                      : "normal",
                  display: "inline-block",
                  overflow: "hidden",
                  lineHeight: "1.1em",
                  height: "1.1em",
                  wordBreak: "break-all",
                  verticalAlign: "middle"
                }}
              >
                {annotations.isEditingLabelName &&
                annotations.labelEditable.category === metadataField &&
                annotations.labelEditable.label === categoryIndex
                  ? null
                  : truncatedString || displayString}
              </span>
            </Tooltip>
            {isUserAnno &&
            annotations.labelEditable.category === metadataField &&
            annotations.isEditingLabelName &&
            annotations.labelEditable.label === categoryIndex ? (
              <form
                onSubmit={e => {
                  e.preventDefault();
                  this.handleEditValue();
                }}
              >
                <InputGroup
                  style={{ position: "relative", top: -1 }}
                  ref={input => {
                    this.editableInput = input;
                  }}
                  small
                  autoFocus
                  onChange={e => {
                    this.setState({ editedLabelText: e.target.value });
                  }}
                  defaultValue={displayString}
                  rightElement={
                    <Button
                      minimal
                      style={{ position: "relative", top: -1 }}
                      type="button"
                      icon="small-tick"
                      data-testclass="submitEdit"
                      data-testid="submitEdit"
                      onClick={this.handleEditValue}
                    />
                  }
                />
              </form>
            ) : null}
            {/*
            CANCEL IT, WITH BUTTON, ESCAPE KEY, CLICK OUT, UNDO?

            <Button
            minimal
            style={{ position: "relative", top: -1 }}
            type="button"
            icon="cross"
            data-testclass="submitEdit"
            data-testid="submitEdit"
            onClick={this.cancelEdit}
          /> */}
          </div>
          <span style={{ flexShrink: 0 }}>
            {colorAccessor && !isColorBy && !annotations.isEditingLabelName ? (
              <Occupancy category={category} {...this.props} />
            ) : null}
          </span>
        </div>
        <span>
          <span
            data-testclass="categorical-value-count"
            data-testid={`categorical-value-count-${metadataField}-${displayString}`}
            style={{
              color:
                displayString === globals.unassignedCategoryLabel
                  ? "#ababab"
                  : "black",
              fontStyle:
                displayString === globals.unassignedCategoryLabel
                  ? "italic"
                  : "auto"
            }}
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
            <span
              onMouseEnter={this.handleMouseExit}
              onMouseLeave={this.handleMouseEnter}
            >
              <Popover
                interactionKind={PopoverInteractionKind.HOVER}
                boundary="window"
                position={Position.RIGHT_TOP}
                content={
                  <Menu>
                    <MenuItem
                      icon="plus"
                      data-testclass="handleAddCurrentSelectionToThisLabel"
                      data-testid={`handleAddCurrentSelectionToThisLabel-${metadataField}`}
                      onClick={this.handleAddCurrentSelectionToThisLabel}
                      text={
                        <span>
                          Re-label currently selected cells as
                          <span
                            style={{
                              fontStyle:
                                displayString ===
                                globals.unassignedCategoryLabel
                                  ? "italic"
                                  : "auto"
                            }}
                          >
                            {` ${displayString}`}
                          </span>
                        </span>
                      }
                      disabled={this.isAddCurrentSelectionDisabled(
                        metadataField,
                        value
                      )}
                    />
                    {displayString !== globals.unassignedCategoryLabel ? (
                      <MenuItem
                        icon="edit"
                        text="Edit this label's name"
                        data-testclass="handleEditValue"
                        data-testid={`handleEditValue-${metadataField}`}
                        onClick={this.activateEditLabelMode}
                      />
                    ) : null}
                    {displayString !== globals.unassignedCategoryLabel ? (
                      <MenuItem
                        icon="delete"
                        intent="danger"
                        data-testclass="handleDeleteValue"
                        data-testid={`handleDeleteValue-${metadataField}`}
                        onClick={this.handleDeleteValue}
                        text={`Delete this label, and reassign all cells to type '${globals.unassignedCategoryLabel}'`}
                      />
                    ) : null}
                  </Menu>
                }
              >
                <Button
                  style={{
                    marginLeft: 0,
                    position: "relative",
                    top: -1,
                    minHeight: 16
                  }}
                  data-testclass="seeActions"
                  data-testid={`seeActions-${metadataField}`}
                  icon="more"
                  small
                  minimal
                />
              </Popover>
            </span>
          ) : null}
        </span>
      </div>
    );
  }
}

export default CategoryValue;
