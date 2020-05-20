import { connect } from "react-redux";
import React from "react";

import {
  Button,
  Menu,
  MenuItem,
  Popover,
  Position,
  Icon,
  PopoverInteractionKind,
  Tooltip,
} from "@blueprintjs/core";
import Occupancy from "./occupancy";
import * as globals from "../../../globals";
import styles from "../categorical.css";
import AnnoDialog from "../annoDialog";
import LabelInput from "../labelInput";

import { AnnotationsHelpers } from "../../../util/stateManager";
import maybeTruncateString from "../../../util/maybeTruncateString";
import { labelPrompt, isLabelErroneous } from "../labelUtil";

/* this is defined outside of the class so we can use it in connect() */
function _currentLabel(ownProps, categoricalSelection) {
  const { metadataField, categoryIndex } = ownProps;
  return String(
    categoricalSelection[metadataField].categoryValues[categoryIndex]
  ).valueOf();
}

@connect((state, ownProps) => {
  const { pointDilation, categoricalSelection } = state;
  const { metadataField } = ownProps;
  const isDilated =
    pointDilation.metadataField === metadataField &&
    pointDilation.categoryField ===
      _currentLabel(ownProps, categoricalSelection);
  return {
    categoricalSelection,
    annotations: state.annotations,
    colorScale: state.colors.scale,
    colorAccessor: state.colors.colorAccessor,
    schema: state.world?.schema,
    world: state.world,
    crossfilter: state.crossfilter,
    ontology: state.ontology,
    isDilated,
  };
})
class CategoryValue extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      editedLabelText: this.currentLabel(),
    };
  }

  componentDidUpdate(prevProps) {
    const { categoricalSelection, metadataField, categoryIndex } = this.props;
    if (
      prevProps.categoricalSelection !== categoricalSelection ||
      prevProps.metadataField !== metadataField ||
      prevProps.categoryIndex !== categoryIndex
    ) {
      // adequately checked to prevent looping
      // eslint-disable-next-line react/no-did-update-set-state
      this.setState({
        editedLabelText: this.currentLabel(),
      });
    }
  }

  handleDeleteValue = () => {
    const { dispatch, metadataField } = this.props;
    const label = this.getLabel();

    dispatch({
      type: "annotation: delete label",
      metadataField,
      label,
    });
  };

  handleAddCurrentSelectionToThisLabel = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    const label = this.getLabel();
    dispatch({
      type: "annotation: label current cell selection",
      metadataField,
      categoryIndex,
      label,
    });
  };

  handleEditValue = (e) => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    const { editedLabelText } = this.state;
    const label = this.getLabel();
    this.cancelEditMode();
    dispatch({
      type: "annotation: label edited",
      editedLabel: editedLabelText,
      metadataField,
      categoryIndex,
      label,
    });
    e.preventDefault();
  };

  handleCreateArbitraryLabel = (txt) => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    const label = this.getLabel();
    this.cancelEditMode();
    dispatch({
      type: "annotation: label edited",
      metadataField,
      editedLabel: txt,
      categoryIndex,
      label,
    });
  };

  labelNameError = (name) => {
    const { metadataField, ontology, schema } = this.props;
    if (name === this.currentLabel()) return false;
    return isLabelErroneous(name, metadataField, ontology, schema);
  };

  instruction = (label) => {
    return labelPrompt(this.labelNameError(label), "New, unique label", ":");
  };

  activateEditLabelMode = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "annotation: activate edit label mode",
      metadataField,
      categoryIndex,
    });
  };

  cancelEditMode = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    this.setState({
      editedLabelText: this.currentLabel(),
    });
    dispatch({
      type: "annotation: cancel edit label mode",
      metadataField,
      categoryIndex,
    });
  };

  toggleOff = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "categorical metadata filter deselect",
      metadataField,
      categoryIndex,
    });
  };

  shouldComponentUpdate = (nextProps, nextState) => {
    /*
    Checks to see if at least one of the following changed:
    * world state
    * the color accessor (what is currently being colored by)
    * if this catagorical value's selection status has changed
    * the crossfilter (ie, global selection state)

    If and only if true, update the component
    */
    const { props, state } = this;
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
    const crossfilterChange =
      props.isUserAnno && props.crossfilter !== nextProps.crossfilter;
    const editingLabel = state.editedLabelText !== nextState.editedLabelText;
    const dilationChange = props.isDilated !== nextProps.isDilated;

    return (
      valueSelectionChange ||
      worldChange ||
      colorAccessorChange ||
      annotationsChange ||
      crossfilterChange ||
      editingLabel ||
      dilationChange
    );
  };

  toggleOn = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "categorical metadata filter select",
      metadataField,
      categoryIndex,
    });
  };

  handleMouseEnter = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "category value mouse hover start",
      metadataField,
      categoryIndex,
    });
  };

  handleMouseExit = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    dispatch({
      type: "category value mouse hover end",
      metadataField,
      categoryIndex,
    });
  };

  handleTextChange = (text) => {
    this.setState({ editedLabelText: text });
  };

  handleChoice = (e) => {
    /* Blueprint Suggest format */
    this.setState({ editedLabelText: e.target });
  };

  getLabel = () => {
    const { metadataField, categoryIndex, categoricalSelection } = this.props;
    const category = categoricalSelection[metadataField];
    const label = category.categoryValues[categoryIndex];

    return label;
  };

  currentLabel() {
    const { categoricalSelection } = this.props;
    return _currentLabel(this.props, categoricalSelection);
  }

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
      ontology,
      // flippedProps is potentially brittle, their docs want {...flippedProps} on our div,
      // our lint doesn't like jsx spread, we are version pinned to prevent api change on their part
      flippedProps,
      isDilated,
      world,
    } = this.props;
    const ontologyEnabled = ontology?.enabled ?? false;

    const { editedLabelText } = this.state;

    if (!categoricalSelection) return null;

    const category = categoricalSelection[metadataField];
    const selected = category.categoryValueSelected[categoryIndex];
    const count = category.categoryValueCounts[categoryIndex];
    const value = category.categoryValues[categoryIndex];
    const displayString = this.currentLabel();

    /* this is the color scale, so add swatches below */
    const isColorBy = metadataField === colorAccessor;
    let categories = null;

    if (isColorBy && schema) {
      categories = schema.annotations.obsByName[colorAccessor]?.categories;
    }

    const truncatedString = maybeTruncateString(
      displayString,
      colorAccessor && !isColorBy
        ? globals.categoryLabelDisplayStringShortLength
        : globals.categoryLabelDisplayStringLongLength
    );

    const editModeActive =
      isUserAnno &&
      annotations.labelEditable.category === metadataField &&
      annotations.isEditingLabelName &&
      annotations.labelEditable.label === categoryIndex;

    const valueToggleLabel = `value-toggle-checkbox-${displayString}`;

    return (
      <div
        key={i}
        data-flip-config={flippedProps["data-flip-config"]}
        data-flip-id={flippedProps["data-flip-id"]}
        data-portal-key={flippedProps["data-portal-key"]}
        className={
          /* This code is to change the styles on centroid label hover is causing over-rendering */
          `${styles.value}${isDilated ? ` ${styles.hover}` : ""}`
        }
        data-testclass="categorical-row"
        style={{
          padding: "4px 0px 4px 7px",
          display: "flex",
          alignItems: "baseline",
          justifyContent: "space-between",
          marginBottom: "2px",
          borderRadius: "2px",
        }}
        onMouseEnter={this.handleMouseEnter}
        onMouseLeave={this.handleMouseExit}
      >
        <div
          style={{
            margin: 0,
            padding: 0,
            userSelect: "none",
            width: globals.leftSidebarWidth - 145,
            display: "flex",
            justifyContent: "space-between",
          }}
        >
          <div style={{ display: "flex", alignItems: "baseline" }}>
            <label
              htmlFor={valueToggleLabel}
              className="bp3-control bp3-checkbox"
              style={{ margin: 0 }}
            >
              <input
                id={valueToggleLabel}
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
              position={Position.LEFT}
              usePortal
              modifiers={{
                preventOverflow: { enabled: false },
                hide: { enabled: false },
              }}
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
                  verticalAlign: "middle",
                }}
              >
                {truncatedString || displayString}
              </span>
            </Tooltip>
            {editModeActive ? (
              <div>
                <AnnoDialog
                  isActive={editModeActive}
                  inputProps={{
                    "data-testid": `${metadataField}:edit-label-name-dialog`,
                  }}
                  primaryButtonProps={{
                    "data-testid": `${metadataField}:${displayString}:submit-label-edit`,
                  }}
                  title="Edit label"
                  instruction={this.instruction(editedLabelText)}
                  cancelTooltipContent="Close this dialog without editing label text."
                  primaryButtonText="Change label text"
                  text={editedLabelText}
                  categoryToDuplicate={null}
                  validationError={this.labelNameError(editedLabelText)}
                  handleSubmit={this.handleEditValue}
                  handleCancel={this.cancelEditMode}
                  annoInput={
                    <LabelInput
                      label={editedLabelText}
                      labelSuggestions={ontologyEnabled ? ontology.terms : null}
                      onChange={this.handleTextChange}
                      onSelect={this.handleTextChange}
                      inputProps={{
                        "data-testid": `${metadataField}:${displayString}:edit-label-name`,
                        leftIcon: "tag",
                        intent: "none",
                        autoFocus: true,
                      }}
                    />
                  }
                  annoSelect={null}
                />
              </div>
            ) : null}
          </div>
          <span style={{ flexShrink: 0 }}>
            {colorAccessor && !isColorBy && !annotations.isEditingLabelName ? (
              <Occupancy
                categoryValue={value}
                colorAccessor={colorAccessor}
                metadataField={metadataField}
                world={world}
                colorScale={colorScale}
                colorByIsCategorical={!!categoricalSelection[colorAccessor]}
              />
            ) : null}
          </span>
        </div>
        <div>
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
                    : "auto",
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
                    : "inherit",
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
                        data-testid={`${metadataField}:${displayString}:add-current-selection-to-this-label`}
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
                                    : "auto",
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
                          data-testid={`${metadataField}:${displayString}:edit-label`}
                          onClick={this.activateEditLabelMode}
                          disabled={annotations.isEditingLabelName}
                        />
                      ) : null}
                      {displayString !== globals.unassignedCategoryLabel ? (
                        <MenuItem
                          icon="delete"
                          intent="danger"
                          data-testclass="handleDeleteValue"
                          data-testid={`${metadataField}:${displayString}:delete-label`}
                          onClick={this.handleDeleteValue}
                          text={`Delete this label, and reassign all cells to type '${globals.unassignedCategoryLabel}'`}
                        />
                      ) : null}
                    </Menu>
                  }
                >
                  <Button
                    style={{
                      marginLeft: 2,
                      position: "relative",
                      top: -1,
                      minHeight: 16,
                    }}
                    data-testclass="seeActions"
                    data-testid={`${metadataField}:${displayString}:see-actions`}
                    icon={<Icon icon="more" iconSize={10} />}
                    small
                    minimal
                  />
                </Popover>
              </span>
            ) : null}
          </span>
        </div>
      </div>
    );
  }
}

export default CategoryValue;
