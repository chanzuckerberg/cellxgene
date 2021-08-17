import React, { useRef, useEffect } from "react";
import { connect, shallowEqual } from "react-redux";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import {
  AnchorButton,
  Button,
  Classes,
  Position,
  Tooltip,
} from "@blueprintjs/core";
import { Flipper, Flipped } from "react-flip-toolkit";
import Async from "react-async";
import memoize from "memoize-one";

import Value from "../value";
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
import { Dataframe } from "../../../util/dataframe";

const LABEL_WIDTH = globals.leftSidebarWidth - 100;
const ANNO_BUTTON_WIDTH = 50;
const LABEL_WIDTH_ANNO = LABEL_WIDTH - ANNO_BUTTON_WIDTH;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state, ownProps) => {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const schema = (state as any).annoMatrix?.schema;
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type '{... Remove this comment to see the full error message
  const { metadataField } = ownProps;
  const isUserAnno = schema?.annotations?.obsByName[metadataField]?.writable;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const categoricalSelection = (state as any).categoricalSelection?.[
    metadataField
  ];
  return {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    colors: (state as any).colors,
    categoricalSelection,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    annotations: (state as any).annotations,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    annoMatrix: (state as any).annoMatrix,
    schema,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    crossfilter: (state as any).obsCrossfilter,
    isUserAnno,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    genesets: (state as any).genesets.genesets,
  };
})
class Category extends React.PureComponent {
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  static getSelectionState(
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    categoricalSelection: any,
    // @ts-expect-error ts-migrate(6133) FIXME: 'metadataField' is declared but its value is never... Remove this comment to see the full error message
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    metadataField: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    categorySummary: any
  ) {
    // total number of categories in this dimension
    const totalCatCount = categorySummary.numCategoryValues;
    // number of selected options in this category
    const selectedCatCount = categorySummary.categoryValues.reduce(
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (res: any, label: any) =>
        categoricalSelection.get(label) ?? true ? res + 1 : res,
      0
    );
    return selectedCatCount === totalCatCount
      ? "all"
      : selectedCatCount === 0
      ? "none"
      : "some";
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  static watchAsync(props: any, prevProps: any) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  createCategorySummaryFromDfCol = memoize(createCategorySummaryFromDfCol);

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  getSelectionState(categorySummary: any) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoricalSelection' does not exist on ... Remove this comment to see the full error message
    const { categoricalSelection, metadataField } = this.props;
    return Category.getSelectionState(
      categoricalSelection,
      metadataField,
      categorySummary
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleColorChange = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "color by categorical metadata",
      colorAccessor: metadataField,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleCategoryClick = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annotations' does not exist on type 'Rea... Remove this comment to see the full error message
    const { annotations, metadataField, onExpansionChange } = this.props;
    const editingCategory =
      annotations.isEditingCategoryName &&
      annotations.categoryBeingEdited === metadataField;
    if (!editingCategory) {
      onExpansionChange(metadataField);
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleCategoryKeyPress = (e: any) => {
    if (e.key === "Enter") {
      this.handleCategoryClick();
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleToggleAllClick = (categorySummary: any) => {
    const isChecked = this.getSelectionState(categorySummary);
    if (isChecked === "all") {
      this.toggleNone(categorySummary);
    } else {
      this.toggleAll(categorySummary);
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  fetchAsyncProps = async (props: any) => {
    const { annoMatrix, metadataField, colors } = props.watchProps;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'crossfilter' does not exist on type 'Rea... Remove this comment to see the full error message
    const { crossfilter } = this.props;

    const [categoryData, categorySummary, colorData] = await this.fetchData(
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

  async fetchData(
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    annoMatrix: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    metadataField: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    colors: any
  ): Promise<
    [
      Dataframe,
      ReturnType<typeof createCategorySummaryFromDfCol>,
      Dataframe | null
    ]
  > {
    /*
    fetch our data and the color-by data if appropriate, and then build a summary
    of our category and a color table for the color-by annotation.
    */
    const { schema } = annoMatrix;
    const { colorAccessor, colorMode } = colors;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'genesets' does not exist on type 'Readon... Remove this comment to see the full error message
    const { genesets } = this.props;
    let colorDataPromise: Promise<Dataframe | null> = Promise.resolve(null);
    if (colorAccessor) {
      const query = createColorQuery(
        colorMode,
        colorAccessor,
        schema,
        genesets
      );
      if (query) colorDataPromise = annoMatrix.fetch(...query);
    }
    const [categoryData, colorData] = await Promise.all<
      Dataframe,
      Dataframe | null
    >([annoMatrix.fetch("obs", metadataField), colorDataPromise]);

    // our data
    const column = categoryData.icol(0);
    const colSchema = schema.annotations.obsByName[metadataField];
    const categorySummary = this.createCategorySummaryFromDfCol(
      column,
      colSchema
    );
    return [categoryData, categorySummary, colorData];
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types -- - FIXME: disabled temporarily on migrate to TS.
  updateColorTable(colorData: Dataframe|null) {
    // color table, which may be null
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  toggleNone(categorySummary: any) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    dispatch(
      actions.selectCategoricalAllMetadataAction(
        "categorical metadata filter none of these",
        metadataField,
        categorySummary.allCategoryValues,
        false
      )
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  toggleAll(categorySummary: any) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    dispatch(
      actions.selectCategoricalAllMetadataAction(
        "categorical metadata filter all of these",
        metadataField,
        categorySummary.allCategoryValues,
        true
      )
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
      metadataField,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isExpanded' does not exist on type 'Read... Remove this comment to see the full error message
      isExpanded,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoricalSelection' does not exist on ... Remove this comment to see the full error message
      categoricalSelection,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'crossfilter' does not exist on type 'Rea... Remove this comment to see the full error message
      crossfilter,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colors' does not exist on type 'Readonly... Remove this comment to see the full error message
      colors,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
      annoMatrix,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isUserAnno' does not exist on type 'Read... Remove this comment to see the full error message
      isUserAnno,
    } = this.props;

    const checkboxID = `category-select-${metadataField}`;

    return (
      <CategoryCrossfilterContext.Provider value={crossfilter}>
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
          <Async.Pending initial>
            <StillLoading
              metadataField={metadataField}
              checkboxID={checkboxID}
            />
          </Async.Pending>
          <Async.Rejected>
            {(error) => (
              <ErrorLoading metadataField={metadataField} error={error} />
            )}
          </Async.Rejected>
          <Async.Fulfilled persist>
            {(asyncProps) => {
              const {
                // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorAccessor' does not exist on type 'u... Remove this comment to see the full error message
                colorAccessor,
                // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorTable' does not exist on type 'unkn... Remove this comment to see the full error message
                colorTable,
                // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorData' does not exist on type 'unkno... Remove this comment to see the full error message
                colorData,
                // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryData' does not exist on type 'un... Remove this comment to see the full error message
                categoryData,
                // @ts-expect-error ts-migrate(2339) FIXME: Property 'categorySummary' does not exist on type ... Remove this comment to see the full error message
                categorySummary,
                // @ts-expect-error ts-migrate(2339) FIXME: Property 'isColorAccessor' does not exist on type ... Remove this comment to see the full error message
                isColorAccessor,
                // @ts-expect-error ts-migrate(2339) FIXME: Property 'handleCategoryToggleAllClick' does not e... Remove this comment to see the full error message
                handleCategoryToggleAllClick,
              } = asyncProps;
              const selectionState = this.getSelectionState(categorySummary);
              return (
                <CategoryRender
                  // @ts-expect-error ts-migrate(2322) FIXME: Type '{ metadataField: any; checkboxID: string; is... Remove this comment to see the full error message
                  metadataField={metadataField}
                  checkboxID={checkboxID}
                  isUserAnno={isUserAnno}
                  isExpanded={isExpanded}
                  isColorAccessor={isColorAccessor}
                  selectionState={selectionState}
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
              );
            }}
          </Async.Fulfilled>
        </Async>
      </CategoryCrossfilterContext.Provider>
    );
  }
}

export default Category;

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
const StillLoading = ({ metadataField, checkboxID }: any) => (
  /*
  We are still loading this category, so render a "busy" signal.
  */
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
        <label
          htmlFor={checkboxID}
          className={`${Classes.CONTROL} ${Classes.CHECKBOX}`}
        >
          <input disabled id={checkboxID} checked type="checkbox" />
          <span className={Classes.CONTROL_INDICATOR} />
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
// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
const ErrorLoading = ({ metadataField, error }: any) => {
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
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type '{... Remove this comment to see the full error message
    metadataField,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'checkboxID' does not exist on type '{ ch... Remove this comment to see the full error message
    checkboxID,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isUserAnno' does not exist on type '{ ch... Remove this comment to see the full error message
    isUserAnno,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isColorAccessor' does not exist on type ... Remove this comment to see the full error message
    isColorAccessor,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isExpanded' does not exist on type '{ ch... Remove this comment to see the full error message
    isExpanded,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'selectionState' does not exist on type '... Remove this comment to see the full error message
    selectionState,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onColorChangeClick' does not exist on ty... Remove this comment to see the full error message
    onColorChangeClick,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onCategoryMenuClick' does not exist on t... Remove this comment to see the full error message
    onCategoryMenuClick,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onCategoryMenuKeyPress' does not exist o... Remove this comment to see the full error message
    onCategoryMenuKeyPress,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onCategoryToggleAllClick' does not exist... Remove this comment to see the full error message
    onCategoryToggleAllClick,
  }) => {
    /*
    Render category name and controls (eg, color-by button).
    */
    const checkboxRef = useRef(null);

    useEffect(() => {
      // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
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
          <label
            className={`${Classes.CONTROL} ${Classes.CHECKBOX}`}
            htmlFor={checkboxID}
          >
            <input
              id={checkboxID}
              data-testclass="category-select"
              data-testid={`${metadataField}:category-select`}
              onChange={onCategoryToggleAllClick}
              ref={checkboxRef}
              checked={selectionState === "all"}
              type="checkbox"
            />
            <span className={Classes.CONTROL_INDICATOR} />
          </label>
          <span
            role="menuitem"
            // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'number | ... Remove this comment to see the full error message
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
                // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'number | ... Remove this comment to see the full error message
                tabIndex="-1"
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
        {/* @ts-expect-error ts-migrate(2322) FIXME: Type '{ metadataField: any; }' is not assignable t... Remove this comment to see the full error message */}
        <AnnoDialogEditCategoryName metadataField={metadataField} />
        {/* @ts-expect-error ts-migrate(2322) FIXME: Type '{ metadataField: any; }' is not assignable t... Remove this comment to see the full error message */}
        <AnnoDialogAddLabel metadataField={metadataField} />
        <div>
          <AnnoMenu
            // @ts-expect-error ts-migrate(2322) FIXME: Type '{ metadataField: any; isUserAnno: any; creat... Remove this comment to see the full error message
            metadataField={metadataField}
            isUserAnno={isUserAnno}
            createText="Add a new label to this category"
            editText="Edit this category's name"
            deleteText="Delete this category, all associated labels, and remove all cell assignments"
          />

          <Tooltip
            content="Use as color scale"
            position={Position.LEFT}
            usePortal
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
            modifiers={{
              preventOverflow: { enabled: false },
              hide: { enabled: false },
            }}
          >
            <AnchorButton
              data-testclass="colorby"
              data-testid={`colorby-${metadataField}`}
              onClick={onColorChangeClick}
              active={isColorAccessor}
              intent={isColorAccessor ? "primary" : "none"}
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
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type '{... Remove this comment to see the full error message
    metadataField,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'checkboxID' does not exist on type '{ ch... Remove this comment to see the full error message
    checkboxID,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isUserAnno' does not exist on type '{ ch... Remove this comment to see the full error message
    isUserAnno,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isColorAccessor' does not exist on type ... Remove this comment to see the full error message
    isColorAccessor,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isExpanded' does not exist on type '{ ch... Remove this comment to see the full error message
    isExpanded,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'selectionState' does not exist on type '... Remove this comment to see the full error message
    selectionState,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryData' does not exist on type '{ ... Remove this comment to see the full error message
    categoryData,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'categorySummary' does not exist on type ... Remove this comment to see the full error message
    categorySummary,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorAccessor' does not exist on type '{... Remove this comment to see the full error message
    colorAccessor,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorData' does not exist on type '{ chi... Remove this comment to see the full error message
    colorData,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorTable' does not exist on type '{ ch... Remove this comment to see the full error message
    colorTable,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onColorChangeClick' does not exist on ty... Remove this comment to see the full error message
    onColorChangeClick,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onCategoryMenuClick' does not exist on t... Remove this comment to see the full error message
    onCategoryMenuClick,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onCategoryMenuKeyPress' does not exist o... Remove this comment to see the full error message
    onCategoryMenuKeyPress,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onCategoryToggleAllClick' does not exist... Remove this comment to see the full error message
    onCategoryToggleAllClick,
  }) => {
    /*
    Render the core of the category, including checkboxes, controls, etc.
    */
    const { numCategoryValues } = categorySummary;
    const isSingularValue = !isUserAnno && numCategoryValues === 1;

    if (isSingularValue) {
      /*
      Entire category has a single value, special case.
      */
      return null;
    }

    /*
    Otherwise, our normal multi-layout layout
    */
    return (
      <div
        style={{
          maxWidth: globals.maxControlsWidth,
        }}
        data-testclass="category"
        data-testid={`category-${metadataField}`}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "baseline",
          }}
        >
          <CategoryHeader
            // @ts-expect-error ts-migrate(2322) FIXME: Type '{ metadataField: any; checkboxID: any; isUse... Remove this comment to see the full error message
            metadataField={metadataField}
            checkboxID={checkboxID}
            isUserAnno={isUserAnno}
            isExpanded={isExpanded}
            isColorAccessor={isColorAccessor}
            selectionState={selectionState}
            onColorChangeClick={onColorChangeClick}
            onCategoryToggleAllClick={onCategoryToggleAllClick}
            onCategoryMenuClick={onCategoryMenuClick}
            onCategoryMenuKeyPress={onCategoryMenuKeyPress}
          />
        </div>
        <div style={{ marginLeft: 26 }}>
          {
            /* values*/
            isExpanded ? (
              <CategoryValueList
                // @ts-expect-error ts-migrate(2322) FIXME: Type '{ isUserAnno: any; metadataField: any; categ... Remove this comment to see the full error message
                isUserAnno={isUserAnno}
                metadataField={metadataField}
                categoryData={categoryData}
                categorySummary={categorySummary}
                colorAccessor={colorAccessor}
                colorData={colorData}
                colorTable={colorTable}
              />
            ) : null
          }
        </div>
      </div>
    );
  }
);

const CategoryValueList = React.memo(
  ({
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isUserAnno' does not exist on type '{ ch... Remove this comment to see the full error message
    isUserAnno,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type '{... Remove this comment to see the full error message
    metadataField,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryData' does not exist on type '{ ... Remove this comment to see the full error message
    categoryData,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'categorySummary' does not exist on type ... Remove this comment to see the full error message
    categorySummary,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorAccessor' does not exist on type '{... Remove this comment to see the full error message
    colorAccessor,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorData' does not exist on type '{ chi... Remove this comment to see the full error message
    colorData,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorTable' does not exist on type '{ ch... Remove this comment to see the full error message
    colorTable,
  }) => {
    const tuples = [...categorySummary.categoryValueIndices];

    /*
    Render the value list.  If this is a user annotation, we use a flipper
    animation, if read-only, we don't bother and save a few bits of perf.
    */
    if (!isUserAnno) {
      return (
        <>
          {tuples.map(([value, index]) => (
            <Value
              key={value}
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ key: any; isUserAnno: any; metadataField: ... Remove this comment to see the full error message
              isUserAnno={isUserAnno}
              metadataField={metadataField}
              categoryIndex={index}
              categoryData={categoryData}
              categorySummary={categorySummary}
              colorAccessor={colorAccessor}
              colorData={colorData}
              colorTable={colorTable}
            />
          ))}
        </>
      );
    }

    /* User annotation */
    const flipKey = tuples.map((t) => t[0]).join("");
    return (
      <Flipper flipKey={flipKey}>
        {tuples.map(([value, index]) => (
          <Flipped key={value} flipId={value}>
            <Value
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ isUserAnno: any; metadataField: any; categ... Remove this comment to see the full error message
              isUserAnno={isUserAnno}
              metadataField={metadataField}
              categoryIndex={index}
              categoryData={categoryData}
              categorySummary={categorySummary}
              colorAccessor={colorAccessor}
              colorData={colorData}
              colorTable={colorTable}
            />
          </Flipped>
        ))}
      </Flipper>
    );
  }
);
