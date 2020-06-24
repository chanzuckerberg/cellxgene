import { connect } from "react-redux";
import React from "react";
import * as d3 from "d3";

import {
  Button,
  Menu,
  MenuItem,
  Popover,
  Position,
  Icon,
  PopoverInteractionKind,
} from "@blueprintjs/core";
import * as globals from "../../../globals";
import styles from "../categorical.css";
import AnnoDialog from "../annoDialog";
import LabelInput from "../labelInput";
import Truncate from "../../util/truncate";

import { AnnotationsHelpers } from "../../../util/stateManager";
import { labelPrompt, isLabelErroneous } from "../labelUtil";
import MiniHistogram from "../../miniHistogram";
import MiniStackedBar from "../../miniStackedBar";

const VALUE_HEIGHT = 11;
const CHART_WIDTH = 100;

/* this is defined outside of the class so we can use it in connect() */
function _currentLabelAsString(ownProps) {
  const { categorySummary, categoryIndex } = ownProps;
  // when called as a function, the String() constructor performs type conversion,
  // and returns a primitive string.
  return String(categorySummary.categoryValues[categoryIndex]);
}

@connect((state, ownProps) => {
  const { pointDilation, categoricalSelection } = state;
  const { metadataField } = ownProps;
  const isDilated =
    pointDilation.metadataField === metadataField &&
    pointDilation.categoryField === _currentLabelAsString(ownProps);
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
      editedLabelText: this.currentLabelAsString(),
    };
  }

  componentDidUpdate(prevProps) {
    const {
      categoricalSelection,
      metadataField,
      categoryIndex,
      categorySummary,
    } = this.props;
    if (
      prevProps.categoricalSelection !== categoricalSelection ||
      prevProps.metadataField !== metadataField ||
      prevProps.categoryIndex !== categoryIndex ||
      prevProps.categorySummary !== categorySummary
    ) {
      // eslint-disable-next-line react/no-did-update-set-state --- adequately checked to prevent looping
      this.setState({
        editedLabelText: this.currentLabelAsString(),
      });
    }
  }

  // If coloring by and this isn't the colorAccessor and it isn't being edited
  get shouldRenderStackedBarOrHistogram() {
    const { colorAccessor, isColorBy, annotations } = this.props;

    return colorAccessor && !isColorBy && !annotations.isEditingLabelName;
  }

  getLabel() {
    const { categoryIndex, categorySummary } = this.props;
    const label = categorySummary.categoryValues[categoryIndex];
    return label;
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
    if (name === this.currentLabelAsString()) return false;
    return isLabelErroneous(name, metadataField, ontology, schema);
  };

  instruction = (label) => {
    return labelPrompt(this.labelNameError(label), "New, unique label", ":");
  };

  activateEditLabelMode = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    const label = this.getLabel();
    dispatch({
      type: "annotation: activate edit label mode",
      metadataField,
      categoryIndex,
      label,
    });
  };

  cancelEditMode = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    const label = this.getLabel();
    this.setState({
      editedLabelText: this.currentLabelAsString(),
    });
    dispatch({
      type: "annotation: cancel edit label mode",
      metadataField,
      categoryIndex,
      label,
    });
  };

  toggleOff = () => {
    const {
      dispatch,
      metadataField,
      categoryIndex,
      categorySummary,
    } = this.props;
    const labels = categorySummary.categoryValues;
    const label = labels[categoryIndex];
    dispatch({
      type: "categorical metadata filter deselect",
      metadataField,
      categoryIndex,
      label,
      labels,
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
    const {
      metadataField,
      categoryIndex,
      categoricalSelection,
      categorySummary,
    } = props;
    const {
      categoryIndex: newCategoryIndex,
      categoricalSelection: newCategoricalSelection,
      categorySummary: newCategorySummary,
    } = nextProps;

    const label = categorySummary.categoryValues[categoryIndex];
    const newLabel = newCategorySummary.categoryValues[newCategoryIndex];
    const labelChanged = label !== newLabel;
    const valueSelectionChange =
      categoricalSelection[metadataField].get(label) !==
      newCategoricalSelection[metadataField].get(newLabel);

    const worldChange = props.world !== nextProps.world;
    const colorAccessorChange = props.colorAccessor !== nextProps.colorAccessor;
    const annotationsChange = props.annotations !== nextProps.annotations;
    const crossfilterChange =
      props.isUserAnno && props.crossfilter !== nextProps.crossfilter;
    const editingLabel = state.editedLabelText !== nextState.editedLabelText;
    const dilationChange = props.isDilated !== nextProps.isDilated;

    const count = categorySummary.categoryValueCounts[categoryIndex];
    const newCount = newCategorySummary.categoryValueCounts[newCategoryIndex];
    const countChanged = count !== newCount;

    return (
      labelChanged ||
      valueSelectionChange ||
      worldChange ||
      colorAccessorChange ||
      annotationsChange ||
      crossfilterChange ||
      editingLabel ||
      dilationChange ||
      countChanged
    );
  };

  toggleOn = () => {
    const {
      dispatch,
      metadataField,
      categoryIndex,
      categorySummary,
    } = this.props;
    const labels = categorySummary.categoryValues;
    const label = labels[categoryIndex];
    dispatch({
      type: "categorical metadata filter select",
      metadataField,
      categoryIndex,
      label,
      labels,
    });
  };

  handleMouseEnter = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    const label = this.getLabel();
    dispatch({
      type: "category value mouse hover start",
      metadataField,
      categoryIndex,
      label,
    });
  };

  handleMouseExit = () => {
    const { dispatch, metadataField, categoryIndex } = this.props;
    const label = this.getLabel();
    dispatch({
      type: "category value mouse hover end",
      metadataField,
      categoryIndex,
      label,
    });
  };

  handleTextChange = (text) => {
    this.setState({ editedLabelText: text });
  };

  handleChoice = (e) => {
    /* Blueprint Suggest format */
    this.setState({ editedLabelText: e.target });
  };

  createHistogramBins = (
    world,
    metadataField,
    colorAccessor,
    value,
    width,
    height
  ) => {
    /*
      Knowing that colorScale is based off continuous data,
      createHistogramBins fetches the continuous data in relation to the cells relevant to the category value.
      It then separates that data into 50 bins for drawing the mini-histogram
    */
    const groupBy = world.obsAnnotations.col(metadataField);

    const col =
      world.obsAnnotations.col(colorAccessor) ||
      world.varData.col(colorAccessor);

    const range = col.summarize();

    const histogramMap = col.histogram(
      50,
      [range.min, range.max],
      groupBy
    ); /* Because the signature changes we really need different names for histogram to differentiate signatures  */

    const bins = histogramMap.has(value)
      ? histogramMap.get(value)
      : new Array(50).fill(0);

    const xScale = d3.scaleLinear().domain([0, bins.length]).range([0, width]);

    const largestBin = Math.max(...bins);

    const yScale = d3.scaleLinear().domain([0, largestBin]).range([0, height]);

    return {
      xScale,
      yScale,
      bins,
    };
  };

  createStackedGraphBins = (
    world,
    metadataField,
    colorAccessor,
    categoryValue,
    width
  ) => {
    /*
      Knowing that the color scale is based off of categorical data,
      createOccupancyStack obtains a map showing the number if cells per colored value
      Using the colorScale a stack of colored bars is drawn representing the map
     */
    const { schema } = world;

    const groupBy = world.obsAnnotations.col(metadataField);
    const occupancyMap = world.obsAnnotations
      .col(colorAccessor)
      .histogramCategorical(groupBy);

    const occupancy = occupancyMap.get(categoryValue);

    if (occupancy && occupancy.size > 0) {
      // not all categories have occupancy, so occupancy may be undefined.
      const scale = d3
        .scaleLinear()
        /* get all the keys d[1] as an array, then find the sum */
        .domain([0, d3.sum(Array.from(occupancy.values()))])
        .range([0, width]);
      const categories =
        schema.annotations.obsByName[colorAccessor]?.categories;

      const dfColumn = world.obsAnnotations.col(colorAccessor);
      const categoryValues = dfColumn.summarizeCategorical().categories;

      return {
        domainValues: categoryValues,
        scale,
        domain: categories,
        occupancy,
      };
    }
    return null;
  };

  currentLabelAsString() {
    return _currentLabelAsString(this.props);
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

  renderMiniStackedBar = (categoryValue) => {
    const {
      categoricalSelection,
      colorAccessor,
      colorScale,
      metadataField,
      world,
    } = this.props;

    const isColorBy = metadataField === colorAccessor;

    if (
      !this.shouldRenderStackedBarOrHistogram ||
      !categoricalSelection[colorAccessor] ||
      isColorBy
    ) {
      return null;
    }

    const { domainValues, scale, domain, occupancy } =
      this.createStackedGraphBins(
        world,
        metadataField,
        colorAccessor,
        categoryValue,
        CHART_WIDTH
      ) ?? {};

    if (!domainValues || !scale || !domain || !occupancy) {
      return null;
    }

    return (
      <MiniStackedBar
        /* eslint-disable react/jsx-props-no-spreading -- Disable unneeded on next release of eslint-config-airbnb */
        {...{
          colorScale,
          domainValues,
          scale,
          domain,
          occupancy,
        }}
        /* eslint-enable react/jsx-props-no-spreading -- enable */
        height={VALUE_HEIGHT}
        width={CHART_WIDTH}
      />
    );
  };

  renderMiniHistogram = (categoryValue) => {
    const {
      categoricalSelection,
      colorAccessor,
      colorScale,
      world,
      metadataField,
    } = this.props;

    if (
      !this.shouldRenderStackedBarOrHistogram ||
      categoricalSelection[colorAccessor]
    ) {
      return null;
    }

    const { xScale, yScale, bins } =
      this.createHistogramBins(
        world,
        metadataField,
        colorAccessor,
        categoryValue,
        CHART_WIDTH,
        VALUE_HEIGHT
      ) ?? {};

    return (
      <MiniHistogram
        /* eslint-disable react/jsx-props-no-spreading -- Disable unneeded on next release of eslint-config-airbnb */
        {...{
          colorScale,
          xScale,
          yScale,
          bins,
        }}
        /* eslint-enable react/jsx-props-no-spreading -- enable */
        obsOrVarContinuousFieldDisplayName={colorAccessor}
        domainLabel={categoryValue}
        height={VALUE_HEIGHT}
        width={CHART_WIDTH}
      />
    );
  };

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
      categorySummary,
    } = this.props;
    const ontologyEnabled = ontology?.enabled ?? false;

    const { editedLabelText } = this.state;

    if (!categoricalSelection) return null;

    const category = categoricalSelection[metadataField];
    const selected = category.get(this.getLabel()) ?? true;
    const count = categorySummary.categoryValueCounts[categoryIndex];
    const value = categorySummary.categoryValues[categoryIndex];
    const displayString = this.currentLabelAsString();

    /* this is the color scale, so add swatches below */
    const isColorBy = metadataField === colorAccessor;
    let categories = null;

    if (isColorBy && schema) {
      categories = schema.annotations.obsByName[colorAccessor]?.categories;
    }

    const editModeActive =
      isUserAnno &&
      annotations.labelEditable.category === metadataField &&
      annotations.isEditingLabelName &&
      annotations.labelEditable.label === categoryIndex;

    const valueToggleLabel = `value-toggle-checkbox-${displayString}`;

    const LEFT_MARGIN = 60;
    const CHECKBOX = 26;
    const CELL_NUMBER = 50;
    const ANNO_MENU = 26;
    const LABEL_MARGIN = 16;
    const CHART_MARGIN = 24;

    const otherElementsWidth =
      LEFT_MARGIN +
      CHECKBOX +
      CELL_NUMBER +
      LABEL_MARGIN +
      (isUserAnno ? ANNO_MENU : 0);

    const labelWidth =
      colorAccessor && !isColorBy
        ? globals.leftSidebarWidth -
          otherElementsWidth -
          CHART_WIDTH -
          CHART_MARGIN
        : globals.leftSidebarWidth - otherElementsWidth;

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
            <Truncate>
              <span
                data-testid={`categorical-value-${metadataField}-${displayString}`}
                data-testclass="categorical-value"
                style={{
                  width: labelWidth,
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
                  verticalAlign: "middle",
                  marginRight: LABEL_MARGIN,
                }}
              >
                {displayString}
              </span>
            </Truncate>
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
            {this.renderMiniStackedBar(value)}
            {this.renderMiniHistogram(value)}
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
                width: VALUE_HEIGHT,
                height: VALUE_HEIGHT,
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
