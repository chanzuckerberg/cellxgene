import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import { AnchorButton, Button, Tooltip } from "@blueprintjs/core";
import CategoryFlipperLayout from "./categoryFlipperLayout";
import AnnoMenu from "./annoMenuCategory";
import AnnoDialogEditCategoryName from "./annoDialogEditCategoryName";
import AnnoDialogAddLabel from "./annoDialogAddLabel";
import Truncate from "../../util/truncate";

import * as globals from "../../../globals";

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
    schema: state.world?.schema,
  };
})
class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isChecked: true,
    };
  }

  componentDidUpdate(prevProps) {
    const { categoricalSelection, metadataField } = this.props;
    const cat = categoricalSelection?.[metadataField];
    if (
      categoricalSelection !== prevProps.categoricalSelection &&
      !!cat &&
      !!this.checkbox
    ) {
      const categoryCount = {
        // total number of categories in this dimension
        totalCatCount: cat.numCategoryValues,
        // number of selected options in this category
        selectedCatCount: _.reduce(
          cat.categoryValueSelected,
          (res, cond) => (cond ? res + 1 : res),
          0
        ),
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
        this.setState({ isChecked: false }); // eslint-disable-line react/no-did-update-set-state
      }
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

  toggleNone() {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "categorical metadata filter none of these",
      metadataField,
    });
    this.setState({ isChecked: false });
  }

  toggleAll() {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "categorical metadata filter all of these",
      metadataField,
    });
    this.setState({ isChecked: true });
  }

  handleToggleAllClick() {
    const { isChecked } = this.state;
    // || this.checkbox.indeterminate === false
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
    const { isChecked } = this.state;
    const {
      metadataField,
      categoricalSelection,
      isColorAccessor,
      isExpanded,
      schema,
    } = this.props;

    const isStillLoading = !(categoricalSelection?.[metadataField] ?? false);
    if (isStillLoading) {
      return this.renderIsStillLoading();
    }

    const checkboxID = `category-select-${metadataField}`;

    const isUserAnno =
      schema?.annotations?.obsByName[metadataField]?.writable ?? false;
    const isTruncated = _.get(
      categoricalSelection,
      [metadataField, "isTruncated"],
      false
    );

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
