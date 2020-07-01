import React, { useRef, useEffect } from "react";
import { connect } from "react-redux";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import { AnchorButton, Button, Tooltip } from "@blueprintjs/core";
import Async from "react-async";
import memoize from "memoize-one";

import CategoryFlipperLayout from "./categoryFlipperLayout";
import AnnoMenu from "./annoMenuCategory";
import AnnoDialogEditCategoryName from "./annoDialogEditCategoryName";
import AnnoDialogAddLabel from "./annoDialogAddLabel";
import Truncate from "../../util/truncate";
import { CategoryCrossfilterContext } from "../categoryContext";

import * as globals from "../../../globals";
import { createCategorySummaryFromDfCol } from "../../../util/stateManager/controlsHelpers";
import {
  createColorTable,
  createColorQuery,
} from "../../../util/stateManager/colorHelpers";
import actions from "../../../actions";
import shallowEqual from "../../../util/shallowEqual";

const LABEL_WIDTH = globals.leftSidebarWidth - 100;
const ANNO_BUTTON_WIDTH = 50;
const LABEL_WIDTH_ANNO = LABEL_WIDTH - ANNO_BUTTON_WIDTH;

@connect((state) => ({
  colors: state.colors,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  annoMatrix: state.annoMatrix,
  schema: state.annoMatrix?.schema,
  crossfilter: state.obsCrossfilter,
}))
class Category extends React.PureComponent {
  static getSelectionState(
    categoricalSelection,
    metadataField,
    categorySummary
  ) {
    const cat = categoricalSelection?.[metadataField];

    // total number of categories in this dimension
    const totalCatCount = categorySummary.numCategoryValues;
    // number of selected options in this category
    const selectedCatCount = categorySummary.categoryValues.reduce(
      (res, label) => (cat.get(label) ?? true ? res + 1 : res),
      0
    );
    return selectedCatCount === totalCatCount
      ? "all"
      : selectedCatCount === 0
      ? "none"
      : "some";
  }

  static watchAsync(props, prevProps) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  static async fetchData(annoMatrix, metadataField, colors) {
    /*
    fetch our data and the color-by data if appropriate, and then build a summary
    of our category and a color table for the color-by annotation.
    */
    const { schema } = annoMatrix;
    const { colorAccessor, colorMode } = colors;
    let colorDataPromise = Promise.resolve(null);
    if (colorAccessor) {
      const query = createColorQuery(colorMode, colorAccessor, schema);
      if (query) colorDataPromise = annoMatrix.fetch(...query);
    }
    const [categoryData, colorData] = await Promise.all([
      annoMatrix.fetch("obs", metadataField),
      colorDataPromise,
    ]);

    // our data
    const column = categoryData.icol(0);
    const colSchema = schema.annotations.obsByName[metadataField];
    const categorySummary = createCategorySummaryFromDfCol(column, colSchema);
    return [categoryData, categorySummary, colorData];
  }

  getSelectionState = memoize((categorySummary) => {
    const { categoricalSelection, metadataField } = this.props;
    return Category.getSelectionState(
      categoricalSelection,
      metadataField,
      categorySummary
    );
  });

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

  handleCategoryKeyPress = (e) => {
    if (e.key === "Enter") {
      this.handleCategoryClick();
    }
  };

  handleToggleAllClick = (categorySummary) => {
    const isChecked = this.getSelectionState(categorySummary);
    if (isChecked === "all") {
      this.toggleNone(categorySummary);
    } else {
      this.toggleAll(categorySummary);
    }
  };

  fetchAsyncProps = async (props) => {
    const { annoMatrix, metadataField, colors } = props.watchProps;
    const { crossfilter } = this.props;

    const [categoryData, categorySummary, colorData] = await Category.fetchData(
      annoMatrix,
      metadataField,
      colors
    );

    return {
      categoryData,
      categorySummary,
      colorData,
      crossfilter,
      ...this.updateColorTable(colorData),
      handleCategoryToggleAllClick: () =>
        this.handleToggleAllClick(categorySummary),
    };
  };

  updateColorTable(colorData) {
    // color table, which may be null
    const { schema, colors, metadataField } = this.props;
    const { colorAccessor, userColors, colorMode } = colors;
    return {
      isColorAccessor: colorAccessor === metadataField,
      colorAccessor,
      colorMode,
      colorTable: createColorTable(
        colorMode,
        colorAccessor,
        colorData,
        schema,
        userColors
      ),
    };
  }

  toggleNone(categorySummary) {
    const { dispatch, metadataField } = this.props;
    dispatch(
      actions.selectCategoricalAllMetadataAction(
        "categorical metadata filter none of these",
        metadataField,
        categorySummary.categoryValues,
        false
      )
    );
  }

  toggleAll(categorySummary) {
    const { dispatch, metadataField } = this.props;
    dispatch(
      actions.selectCategoricalAllMetadataAction(
        "categorical metadata filter all of these",
        metadataField,
        categorySummary.categoryValues,
        true
      )
    );
  }

  render() {
    const {
      metadataField,
      isExpanded,
      schema,
      categoricalSelection,
      crossfilter,
      colors,
      annoMatrix,
    } = this.props;

    const checkboxID = `category-select-${metadataField}`;
    const isUserAnno = !!schema?.annotations?.obsByName[metadataField]
      ?.writable;

    if (
      !isUserAnno &&
      schema?.annotations?.obsByName[metadataField]?.categories?.length === 1
    ) {
      /*
      Entire category has a single value, special case.
      */
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
      <Async
        watchFn={Category.watchAsync}
        promiseFn={this.fetchAsyncProps}
        watchProps={{
          metadataField,
          annoMatrix,
          categoricalSelection,
          colors,
        }}
      >
        <Async.Pending>
          <StillLoading metadataField={metadataField} checkboxID={checkboxID} />
        </Async.Pending>
        <Async.Rejected>
          {(error) => (
            <ErrorLoading metadataField={metadataField} error={error} />
          )}
        </Async.Rejected>
        <Async.Fulfilled>
          {(asyncProps) => {
            const {
              colorAccessor,
              colorTable,
              colorData,
              categoryData,
              categorySummary,
              isColorAccessor,
              handleCategoryToggleAllClick,
            } = asyncProps;
            return (
              <CategoryCrossfilterContext.Provider value={crossfilter}>
                <CategoryRender
                  metadataField={metadataField}
                  checkboxID={checkboxID}
                  isUserAnno={isUserAnno}
                  isTruncated={!!categorySummary?.isTruncated}
                  isExpanded={isExpanded}
                  isColorAccessor={isColorAccessor}
                  selectionState={this.getSelectionState(categorySummary)}
                  categoryData={categoryData}
                  categorySummary={categorySummary}
                  colorAccessor={colorAccessor}
                  colorData={colorData}
                  colorTable={colorTable}
                  onColorChangeClick={this.handleColorChange}
                  onCategoryToggleAllClick={handleCategoryToggleAllClick}
                  onCategoryMenuClick={this.handleCategoryClick}
                  onCategoryMenuKeyPress={this.handleCategoryKeyPress}
                />
              </CategoryCrossfilterContext.Provider>
            );
          }}
        </Async.Fulfilled>
      </Async>
    );
  }
}

export default Category;

const StillLoading = ({ metadataField, checkboxID }) => {
  /*
  We are still loading this category, so render a "busy" signal.
  */
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
};

const ErrorLoading = ({ metadataField, error }) => {
  console.error(error); // log error to console as it is unexpected.
  return (
    <div style={{ marginBottom: 10, marginTop: 4 }}>
      <span
        style={{
          cursor: "pointer",
          display: "inline-block",
          width: LABEL_WIDTH,
          fontStyle: "italic",
        }}
      >
        {`Failure loading ${metadataField}`}
      </span>
    </div>
  );
};

const CategoryHeader = React.memo(
  ({
    metadataField,
    checkboxID,
    isUserAnno,
    isTruncated,
    isColorAccessor,
    isExpanded,
    selectionState,
    onColorChangeClick,
    onCategoryMenuClick,
    onCategoryMenuKeyPress,
    onCategoryToggleAllClick,
  }) => {
    /*
    Render category name and controls (eg, color-by button).
    */
    const checkboxRef = useRef(null);

    useEffect(() => {
      checkboxRef.current.indeterminate = selectionState === "some";
    }, [checkboxRef.current, selectionState]);

    return (
      <>
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
              onChange={onCategoryToggleAllClick}
              ref={checkboxRef}
              checked={selectionState === "all"}
              type="checkbox"
            />
            <span className="bp3-control-indicator" />
          </label>
          <span
            role="menuitem"
            tabIndex="0"
            data-testclass="category-expand"
            data-testid={`${metadataField}:category-expand`}
            onKeyPress={onCategoryMenuKeyPress}
            style={{
              cursor: "pointer",
            }}
            onClick={onCategoryMenuClick}
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
              onClick={onColorChangeClick}
              active={isColorAccessor}
              intent={isColorAccessor ? "primary" : "none"}
              disabled={isTruncated}
              icon="tint"
            />
          </Tooltip>
        </div>
      </>
    );
  }
);

const CategoryRender = React.memo(
  ({
    metadataField,
    checkboxID,
    isUserAnno,
    isTruncated,
    isColorAccessor,
    isExpanded,
    selectionState,
    categoryData,
    categorySummary,
    colorAccessor,
    colorData,
    colorTable,
    onColorChangeClick,
    onCategoryMenuClick,
    onCategoryMenuKeyPress,
    onCategoryToggleAllClick,
  }) => {
    /*
    Render the core of the category, including checkboxes, controls, etc.
    */

    return (
      <CategoryFlipperLayout
        metadataField={metadataField}
        isExpanded={isExpanded}
        isUserAnno={isUserAnno}
        categoryData={categoryData}
        categorySummary={categorySummary}
        colorAccessor={colorAccessor}
        colorData={colorData}
        colorTable={colorTable}
      >
        <CategoryHeader
          metadataField={metadataField}
          checkboxID={checkboxID}
          isUserAnno={isUserAnno}
          isTruncated={isTruncated}
          isExpanded={isExpanded}
          isColorAccessor={isColorAccessor}
          selectionState={selectionState}
          onColorChangeClick={onColorChangeClick}
          onCategoryToggleAllClick={onCategoryToggleAllClick}
          onCategoryMenuClick={onCategoryMenuClick}
          onCategoryMenuKeyPress={onCategoryMenuKeyPress}
        />
      </CategoryFlipperLayout>
    );
  }
);
