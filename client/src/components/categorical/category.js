import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import { Flipper, Flipped } from "react-flip-toolkit";
import {
  Button,
  Tooltip,
  Menu,
  MenuItem,
  Popover,
  Icon,
  Position,
  PopoverInteractionKind,
  Colors
} from "@blueprintjs/core";
import AnnoDialog from "./annoDialog";
import AnnoInputs from "./annoInputs";

import * as globals from "../../globals";
import Value from "./value";
import { AnnotationsHelpers } from "../../util/stateManager";
import { labelErrorMessage, isLabelErroneous } from "./labelUtil";

@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  ontology: state.ontology,
  ontologyLoading: state.ontology?.loading,
  ontologyEnabled: state.ontology?.enabled
}))
class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isChecked: true,
      isExpanded: false,
      newCategoryText: props.metadataField,
      newLabelText: ""
    };
  }

  componentDidUpdate(prevProps) {
    const { categoricalSelection, metadataField } = this.props;
    if (categoricalSelection !== prevProps.categoricalSelection) {
      const cat = categoricalSelection[metadataField];
      const categoryCount = {
        // total number of categories in this dimension
        totalCatCount: cat.numCategoryValues,
        // number of selected options in this category
        selectedCatCount: _.reduce(
          cat.categoryValueSelected,
          (res, cond) => (cond ? res + 1 : res),
          0
        )
      };
      if (categoryCount.selectedCatCount === categoryCount.totalCatCount) {
        /* everything is on, so not indeterminate */
        this.checkbox.indeterminate = false;
        this.setState({ isChecked: true }); // eslint-disable-line react/no-did-update-set-state
      } else if (categoryCount.selectedCatCount === 0) {
        /* nothing is on, so no */
        this.checkbox.indeterminate = false;
        this.setState({ isChecked: false }); // eslint-disable-line react/no-did-update-set-state
      } else if (categoryCount.selectedCatCount < categoryCount.totalCatCount) {
        /* to be explicit... */
        this.checkbox.indeterminate = true;
        this.setState({ isChecked: false });
      }
    }
  }

  activateAddNewLabelMode = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "annotation: activate add new label mode",
      data: metadataField
    });
  };

  disableAddNewLabelMode = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "annotation: disable add new label mode"
    });
    this.setState({
      newLabelText: ""
    });
  };

  handleAddNewLabelToCategory = () => {
    const { dispatch, metadataField } = this.props;
    const { newLabelText } = this.state;

    dispatch({
      type: "annotation: add new label to category",
      metadataField,
      newLabelText,
      assignSelectedCells: false
    });
    this.setState({ newLabelText: "" });
  };

  addLabelAndAssignCells = () => {
    const { dispatch, metadataField } = this.props;
    const { newLabelText } = this.state;

    dispatch({
      type: "annotation: add new label to category",
      metadataField,
      newLabelText,
      assignSelectedCells: true
    });

    this.setState({ newLabelText: "" });
  };

  handleCreateArbitraryLabel = newLabelTextNotInOntology => {
    const { dispatch, metadataField } = this.props;

    dispatch({
      type: "annotation: add new label to category",
      metadataField,
      newLabelText: newLabelTextNotInOntology,
      assignSelectedCells: false
    });
    this.setState({ newLabelText: "" });
  };

  handleCategoryEditTextChange = txt => {
    this.setState({
      newCategoryText: txt
    });
  };

  activateEditCategoryMode = () => {
    const { dispatch, metadataField } = this.props;

    dispatch({
      type: "annotation: activate category edit mode",
      data: metadataField
    });
  };

  disableEditCategoryMode = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "annotation: disable category edit mode"
    });
  };

  handleEditCategory = () => {
    const { dispatch, metadataField, categoricalSelection } = this.props;
    const { newCategoryText } = this.state;

    const allCategoryNames = _.keys(categoricalSelection);

    if (
      (allCategoryNames.indexOf(newCategoryText) > -1 &&
        newCategoryText !== metadataField) ||
      newCategoryText === ""
    ) {
      return;
    }

    dispatch({
      type: "annotation: category edited",
      metadataField,
      newCategoryText,
      data: newCategoryText
    });
  };

  handleDeleteCategory = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "annotation: delete category",
      metadataField
    });
  };

  handleColorChange = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "color by categorical metadata",
      colorAccessor: metadataField
    });
  };

  labelNameError = name => {
    const { metadataField, ontology, universe } = this.props;
    return isLabelErroneous(name, metadataField, ontology, universe.schema);
  };

  labelNameErrorMessage = name => {
    const { metadataField, ontology, universe } = this.props;
    return labelErrorMessage(name, metadataField, ontology, universe.schema);
  };

  categoryNameErrorMessage = () => {
    const err = this.editedCategoryNameError();
    if (err === false) return null;

    const errorMessageMap = {
      /* map error code to human readable error message */
      "empty-string": "Blank names not allowed",
      duplicate: "Category name must be unique",
      "trim-spaces": "Leading and trailing spaces not allowed",
      "illegal-characters":
        "Only alphanumeric and special characters (-_.) allowed",
      "multi-space-run": "Multiple consecutive spaces not allowed"
    };
    const errorMessage = errorMessageMap[err] ?? "error";
    return (
      <span
        style={{
          display: "block",
          fontStyle: "italic",
          fontSize: 12,
          marginTop: 5,
          color: Colors.ORANGE3
        }}
      >
        {errorMessage}
      </span>
    );
  };

  editedCategoryNameError = () => {
    const { metadataField, categoricalSelection } = this.props;
    const { newCategoryText } = this.state;

    /* check for syntax errors in category name */
    const error = AnnotationsHelpers.annotationNameIsErroneous(newCategoryText);
    if (error) {
      return error;
    }

    /* check for duplicative categories */
    const allCategoryNames = _.keys(categoricalSelection);
    const categoryNameAlreadyExists =
      allCategoryNames.indexOf(newCategoryText) > -1;
    const sameName = newCategoryText === metadataField;
    if (categoryNameAlreadyExists && !sameName) {
      return "duplicate";
    }

    /* otherwise, no error */
    return false;
  };

  /* leaky to have both of these in multiple components */
  handleChoice = e => {
    this.setState({ newLabelText: e.target });
  };

  handleTextChange = text => {
    this.setState({ newLabelText: text });
  };

  toggleAll() {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "categorical metadata filter all of these",
      metadataField
    });
    this.setState({ isChecked: true });
  }

  toggleNone() {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "categorical metadata filter none of these",
      metadataField
    });
    this.setState({ isChecked: false });
  }

  handleToggleAllClick() {
    const { isChecked } = this.state;
    // || this.checkbox.indeterminate === false
    if (isChecked) {
      this.toggleNone();
    } else if (!isChecked) {
      this.toggleAll();
    }
  }

  renderCategoryItems(optTuples) {
    const { metadataField, isUserAnno } = this.props;

    return _.map(optTuples, (tuple, i) => {
      return (
        <Flipped key={tuple[1]} flipId={tuple[1]}>
          {flippedProps => (
            <Value
              isUserAnno={isUserAnno}
              optTuples={optTuples}
              key={tuple[1]}
              metadataField={metadataField}
              categoryIndex={tuple[1]}
              i={i}
              flippedProps={flippedProps}
            />
          )}
        </Flipped>
      );
    });
  }

  render() {
    const { isExpanded, isChecked, newLabelText, newCategoryText } = this.state;
    const {
      metadataField,
      colorAccessor,
      categoricalSelection,
      isUserAnno,
      annotations,
      ontologyEnabled
    } = this.props;
    const { isTruncated } = categoricalSelection[metadataField];
    const cat = categoricalSelection[metadataField];
    const optTuples = [...cat.categoryValueIndices];
    const optTuplesAsKey = _.map(optTuples, t => t[0]).join(""); // animation

    return (
      <div
        style={{
          maxWidth: globals.maxControlsWidth
        }}
        data-testclass="category"
        data-testid={`category-${metadataField}`}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "baseline"
          }}
        >
          <div
            style={{
              display: "flex",
              justifyContent: "flex-start",
              alignItems: "flex-start"
            }}
          >
            <label className="bp3-control bp3-checkbox">
              <input
                data-testclass="category-select"
                data-testid={`category-select-${metadataField}`}
                onChange={this.handleToggleAllClick.bind(this)}
                ref={el => {
                  this.checkbox = el;
                  return el;
                }}
                checked={isChecked}
                type="checkbox"
              />
              <span className="bp3-control-indicator" />
            </label>
            <span
              data-testid={`category-expand-${metadataField}`}
              style={{
                cursor: "pointer",
                display: "inline-block"
              }}
              onClick={() => {
                const editingCategory =
                  annotations.isEditingCategoryName &&
                  annotations.categoryBeingEdited === metadataField;
                if (!editingCategory) {
                  this.setState({ isExpanded: !isExpanded });
                }
              }}
            >
              {isUserAnno ? (
                <Icon style={{ marginRight: 5 }} icon="tag" iconSize={16} />
              ) : null}

              {metadataField}

              <AnnoDialog
                isActive={
                  annotations.isEditingCategoryName &&
                  annotations.categoryBeingEdited === metadataField
                }
                title="Edit category name"
                instruction="New, unique category name:"
                cancelTooltipContent="Close this dialog without editing this category."
                primaryButtonText="Edit category name"
                text={newCategoryText}
                validationError={this.editedCategoryNameError(newCategoryText)}
                errorMessage={this.categoryNameErrorMessage(newCategoryText)}
                handleSubmit={this.handleEditCategory}
                handleCancel={this.disableEditCategoryMode}
                annoInput={
                  <AnnoInputs
                    useSuggest={false}
                    text={newCategoryText}
                    handleTextChange={this.handleCategoryEditTextChange}
                  />
                }
              />

              {isExpanded ? (
                <FaChevronDown
                  data-testclass="category-expand-is-expanded"
                  style={{ fontSize: 10, marginLeft: 5 }}
                />
              ) : (
                <FaChevronRight
                  data-testclass="category-expand-is-not-expanded"
                  style={{ fontSize: 10, marginLeft: 5 }}
                />
              )}
            </span>
          </div>
          <div>
            {isUserAnno ? (
              <>
                <AnnoDialog
                  isActive={
                    annotations.isAddingNewLabel &&
                    annotations.categoryAddingNewLabel === metadataField
                  }
                  title="Add new label to category"
                  instruction="New, unique label name:"
                  cancelTooltipContent="Close this dialog without adding a label."
                  primaryButtonText="Add"
                  secondaryButtonText="Add label & assign currently selected cells"
                  handleSecondaryButtonSubmit={this.addLabelAndAssignCells}
                  text={newLabelText}
                  validationError={this.labelNameError(newLabelText)}
                  errorMessage={this.labelNameErrorMessage(newLabelText)}
                  handleSubmit={this.handleAddNewLabelToCategory}
                  handleCancel={this.disableAddNewLabelMode}
                  annoInput={
                    <AnnoInputs
                      useSuggest={ontologyEnabled}
                      text={newLabelText}
                      handleCreateArbitraryLabel={
                        this.handleCreateArbitraryLabel
                      }
                      handleItemChange={this.handleSuggestActiveItemChange}
                      handleChoice={this.handleChoice}
                      handleTextChange={this.handleTextChange}
                      isTextInvalid={this.labelNameError}
                      isTextInvalidErrorMessage={this.labelNameErrorMessage}
                    />
                  }
                />
                <Popover
                  interactionKind={PopoverInteractionKind.HOVER}
                  boundary="window"
                  position={Position.RIGHT_TOP}
                  content={
                    <Menu>
                      <MenuItem
                        icon="tag"
                        data-testclass="handleAddNewLabelToCategory"
                        data-testid={`handleAddNewLabelToCategory-${metadataField}`}
                        onClick={this.activateAddNewLabelMode}
                        text="Add a new label to this category"
                      />
                      <MenuItem
                        icon="edit"
                        disabled={annotations.isEditingCategoryName}
                        data-testclass="activateEditCategoryMode"
                        data-testid={`activateEditCategoryMode-${metadataField}`}
                        onClick={this.activateEditCategoryMode}
                        text="Edit this category's name"
                      />
                      <MenuItem
                        icon="delete"
                        intent="danger"
                        data-testclass="handleDeleteCategory"
                        data-testid={`handleDeleteCategory-${metadataField}`}
                        onClick={this.handleDeleteCategory}
                        text="Delete this category, all associated labels, and remove all cell assignments"
                      />
                    </Menu>
                  }
                >
                  <Button
                    style={{ marginLeft: 0 }}
                    data-testclass="seeActions"
                    data-testid={`seeActions-${metadataField}`}
                    icon="more"
                    minimal
                  />
                </Popover>
              </>
            ) : null}
            <Tooltip
              content="Use as color scale"
              position="bottom"
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
            >
              <Button
                data-testclass="colorby"
                data-testid={`colorby-${metadataField}`}
                onClick={this.handleColorChange}
                active={colorAccessor === metadataField}
                intent={colorAccessor === metadataField ? "primary" : "none"}
                icon="tint"
              />
            </Tooltip>
          </div>
        </div>
        <div style={{ marginLeft: 26 }}>
          <Flipper spring="veryGentle" flipKey={optTuplesAsKey}>
            {isExpanded ? this.renderCategoryItems(optTuples) : null}
          </Flipper>
        </div>
        <div>
          {isExpanded && isTruncated ? (
            <p style={{ paddingLeft: 15 }}>... truncated list ...</p>
          ) : null}
        </div>
      </div>
    );
  }
}

export default Category;
