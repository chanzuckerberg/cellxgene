import React from "react";
import { connect } from "react-redux";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import { AnchorButton, Button, Tooltip } from "@blueprintjs/core";
import CategoryFlipperLayout from "./categoryFlipperLayout";
import AnnoMenu from "./annoMenuCategory";
import AnnoDialogEditCategoryName from "./annoDialogEditCategoryName";
import AnnoDialogAddLabel from "./annoDialogAddLabel";
import Truncate from "../../util/truncate";

import * as globals from "../../../globals";
import { createCategorySummary as _createCategorySummary } from "../../../util/stateManager/controlsHelpers";

const LABEL_WIDTH = globals.leftSidebarWidth - 100;
const ANNO_BUTTON_WIDTH = 50;
const LABEL_WIDTH_ANNO = LABEL_WIDTH - ANNO_BUTTON_WIDTH;

@connect((state, ownProps) => {
  const { metadataField } = ownProps;
  return {
    isColorAccessor: state.colors.colorAccessor === metadataField,
    categoricalSelection: state.categoricalSelection,
    annotations: state.annotations,
    universe: state.universe,
    world: state.world,
    schema: state.world?.schema,
  };
})
class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isChecked: true,
      categorySummary: this.createCategorySummary(),
    };
  }

  componentDidUpdate(prevProps) {
    const { categoricalSelection, metadataField, world } = this.props;
    let { categorySummary } = this.state;

    if (
      world !== prevProps.world ||
      metadataField !== prevProps.metadataField ||
      !categorySummary
    ) {
      const newCategorySummary = this.createCategorySummary();
      if (categorySummary !== newCategorySummary) {
        categorySummary = newCategorySummary;
        this.setState({ categorySummary }); // eslint-disable-line react/no-did-update-set-state
      }
    }

    const cat = categoricalSelection?.[metadataField];
    if (
      categoricalSelection !== prevProps.categoricalSelection &&
      !!cat &&
      !!this.checkbox
    ) {
      // total number of categories in this dimension
      const totalCatCount = categorySummary.numCategoryValues;
      // number of selected options in this category
      const selectedCatCount = categorySummary.categoryValues.reduce(
        (res, label) => (cat.get(label) ?? true ? res + 1 : res),
        0
      );
      /* eslint-disable react/no-did-update-set-state -- Contained in if statement to prevent infinite looping */
      if (selectedCatCount === totalCatCount) {
        /* everything is on, so not indeterminate */
        this.checkbox.indeterminate = false;
        this.setState({ isChecked: true });
      } else if (selectedCatCount === 0) {
        /* nothing is on, so no */
        this.checkbox.indeterminate = false;
        this.setState({ isChecked: false });
      } else if (selectedCatCount < totalCatCount) {
        /* to be explicit... */
        this.checkbox.indeterminate = true;
        this.setState({ isChecked: false });
      }
      /* eslint-enable react/no-did-update-set-state -- re-enabling*/
    }
  }

  handleColorChange = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "color by categorical metadata",
      colorAccessor: metadataField,
    });
  };

  handleCategoryClick = () => {
    const { annotations, metadataField, onExpansionChange } = this.props;
    const editingCategory =
      annotations.isEditingCategoryName &&
      annotations.categoryBeingEdited === metadataField;
    if (!editingCategory) {
      onExpansionChange(metadataField);
    }
  };

  createCategorySummary() {
    const { world, metadataField } = this.props;
    if (!world || !metadataField || !world.obsAnnotations.hasCol(metadataField))
      return null;
    return _createCategorySummary(world, metadataField);
  }

  toggleNone() {
    const { dispatch, metadataField } = this.props;
    const { categorySummary } = this.state;
    dispatch({
      type: "categorical metadata filter none of these",
      metadataField,
      labels: categorySummary.categoryValues,
    });
    this.setState({ isChecked: false });
  }

  toggleAll() {
    const { dispatch, metadataField } = this.props;
    const { categorySummary } = this.state;
    dispatch({
      type: "categorical metadata filter all of these",
      metadataField,
      labels: categorySummary.categoryValues,
    });
    this.setState({ isChecked: true });
  }

  handleToggleAllClick() {
    const { isChecked } = this.state;
    if (isChecked) {
      this.toggleNone();
    } else {
      this.toggleAll();
    }
  }

  renderIsStillLoading() {
    /*
    We are still loading this category, so render a "busy" signal.
    */
    const { metadataField } = this.props;

    const checkboxID = `category-select-${metadataField}`;

    return (
      <div
        style={{
          maxWidth: globals.maxControlsWidth,
        }}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "baseline",
          }}
        >
          <div
            style={{
              display: "flex",
              justifyContent: "flex-start",
              alignItems: "flex-start",
            }}
          >
            <label htmlFor={checkboxID} className="bp3-control bp3-checkbox">
              <input disabled id={checkboxID} checked type="checkbox" />
              <span className="bp3-control-indicator" />
            </label>
            <Truncate>
              <span
                style={{
                  cursor: "pointer",
                  display: "inline-block",
                  width: LABEL_WIDTH,
                }}
              >
                {metadataField}
              </span>
            </Truncate>
          </div>
          <div>
            <Button minimal loading intent="primary" />
          </div>
        </div>
      </div>
    );
  }

  render() {
    const { isChecked, categorySummary } = this.state;
    const { metadataField, isColorAccessor, isExpanded, schema } = this.props;

    const isStillLoading = !categorySummary;
    if (isStillLoading) {
      return this.renderIsStillLoading();
    }

    const checkboxID = `category-select-${metadataField}`;

    const isUserAnno = !!schema?.annotations?.obsByName[metadataField]
      ?.writable;
    const isTruncated = !!categorySummary?.isTruncated;

    if (
      !isUserAnno &&
      schema?.annotations?.obsByName[metadataField]?.categories?.length === 1
    ) {
      return (
        <div style={{ marginBottom: 10, marginTop: 4 }}>
          <Truncate>
            <span style={{ maxWidth: 150, fontWeight: 700 }}>
              {metadataField}
            </span>
          </Truncate>
          <Truncate>
            <span style={{ maxWidth: 150 }}>
              {`: ${schema.annotations.obsByName[metadataField].categories[0]}`}
            </span>
          </Truncate>
        </div>
      );
    }

    return (
      <CategoryFlipperLayout
        metadataField={metadataField}
        isExpanded={isExpanded}
        isUserAnno={isUserAnno}
        categorySummary={categorySummary}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "flex-start",
            alignItems: "flex-start",
          }}
        >
          <label className="bp3-control bp3-checkbox" htmlFor={checkboxID}>
            <input
              id={checkboxID}
              data-testclass="category-select"
              data-testid={`${metadataField}:category-select`}
              onChange={this.handleToggleAllClick.bind(this)}
              ref={(el) => {
                this.checkbox = el;
                return el;
              }}
              checked={isChecked}
              type="checkbox"
            />
            <span className="bp3-control-indicator" />
          </label>
          <span
            role="menuitem"
            tabIndex="0"
            data-testid={`${metadataField}:category-expand`}
            onKeyPress={(e) => {
              if (e.key === "Enter") {
                this.handleCategoryClick();
              }
            }}
            style={{
              cursor: "pointer",
            }}
            onClick={this.handleCategoryClick}
          >
            <Truncate>
              <span
                style={{
                  maxWidth: isUserAnno ? LABEL_WIDTH_ANNO : LABEL_WIDTH,
                }}
                data-testid={`${metadataField}:category-label`}
              >
                {metadataField}
              </span>
            </Truncate>
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
        {<AnnoDialogEditCategoryName metadataField={metadataField} />}
        {<AnnoDialogAddLabel metadataField={metadataField} />}
        <div>
          <AnnoMenu
            metadataField={metadataField}
            isUserAnno={isUserAnno}
            createText="Add a new label to this category"
            editText="Edit this category's name"
            deleteText="Delete this category, all associated labels, and remove all cell assignments"
          />

          <Tooltip
            content={
              isTruncated
                ? `Coloring by ${metadataField} is disabled, as it exceeds the limit of ${globals.maxCategoricalOptionsToDisplay} labels`
                : "Use as color scale"
            }
            position="bottom"
            usePortal={false}
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
          >
            <AnchorButton
              data-testclass="colorby"
              data-testid={`colorby-${metadataField}`}
              onClick={this.handleColorChange}
              active={isColorAccessor}
              intent={isColorAccessor ? "primary" : "none"}
              disabled={isTruncated}
              icon="tint"
            />
          </Tooltip>
        </div>
      </CategoryFlipperLayout>
    );
  }
}

export default Category;
