import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import { Button, Tooltip, Icon } from "@blueprintjs/core";
import CategoryFlipperLayout from "./categoryFlipperLayout";
import AnnoMenu from "./annoMenuCategory";
import AnnoDialogEditCategoryName from "./annoDialogEditCategoryName";
import AnnoDialogAddLabel from "./annoDialogAddLabel";

import * as globals from "../../globals";

@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe
}))
class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isChecked: true,
      isExpanded: false
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

  handleColorChange = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "color by categorical metadata",
      colorAccessor: metadataField
    });
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

  renderIsStillLoading() {
    /*
    We are still loading this category, so render a "busy" signal.
    */
    const { metadataField } = this.props;
    return (
      <div
        style={{
          maxWidth: globals.maxControlsWidth
        }}
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
              <input disabled checked={true} type="checkbox" />
              <span className="bp3-control-indicator" />
            </label>
            <span
              style={{
                cursor: "pointer",
                display: "inline-block"
              }}
            >
              {metadataField}
            </span>
          </div>
          <div>
            <Button minimal loading intent="primary" />
          </div>
        </div>
      </div>
    );
  }

  render() {
    const { isExpanded, isChecked } = this.state;
    const {
      metadataField,
      categoricalSelection,
      colorAccessor,
      isUserAnno,
      annotations
    } = this.props;

    const isStillLoading = !(categoricalSelection?.[metadataField] ?? false);
    if (isStillLoading) {
      return this.renderIsStillLoading();
    }

    const isTruncated = _.get(categoricalSelection, [metadataField, "isTruncated"], false);

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
            alignItems: "flex-start"
          }}
        >
          <label className="bp3-control bp3-checkbox">
            <input
              data-testclass="category-select"
              data-testid={`${metadataField}:category-select`}
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
            data-testid={`${metadataField}:category-expand`}
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
              ? `Can't color by this category, too many label values (max ${globals.maxCategoricalOptionsToDisplay})`
              : "Use as color scale"
            }
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
          >
            <Button
              data-testclass="colorby"
              data-testid={`colorby-${metadataField}`}
              onClick={() => {
                if (!isTruncated) this.handleColorChange();
              }}
              active={colorAccessor === metadataField || isTruncated}
              intent={colorAccessor === metadataField && !isTruncated ? "primary" : "none"}
              icon="tint"
            />
          </Tooltip>
        </div>
      </CategoryFlipperLayout>
    );
  }
}

export default Category;
