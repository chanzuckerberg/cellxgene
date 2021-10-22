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
import * as d3 from "d3";
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

const LABEL_WIDTH = globals.leftSidebarWidth - 100;
const ANNO_BUTTON_WIDTH = 50;
const LABEL_WIDTH_ANNO = LABEL_WIDTH - ANNO_BUTTON_WIDTH;

@connect((state, ownProps) => {
  const schema = state.annoMatrix?.schema;
  const { metadataField } = ownProps;
  const isUserAnno = schema?.annotations?.obsByName[metadataField]?.writable;
  const categoricalSelection = state.categoricalSelection?.[metadataField];
  return {
    colors: state.colors,
    categoricalSelection,
    annotations: state.annotations,
    annoMatrix: state.annoMatrix,
    schema,
    crossfilter: state.obsCrossfilter,
    isUserAnno,
    layoutChoiceSankey: state.layoutChoice.sankey,
    genesets: state.genesets.genesets,
    sankeySelected: state.sankeySelection.categories?.[metadataField] ?? false,
    displaySankey: state.sankeySelection.displaySankey
  };
})
class Category extends React.PureComponent {
  constructor(props){
    super(props);
    this.state = {
      sortDirection: null
    }
  }
  static getSelectionState(
    categoricalSelection,
    metadataField,
    categorySummary
  ) {
    // total number of categories in this dimension
    const totalCatCount = categorySummary.numCategoryValues;
    // number of selected options in this category
    const selectedCatCount = categorySummary.categoryValues.reduce(
      (res, label) => (categoricalSelection.get(label) ?? true ? res + 1 : res),
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

  createCategorySummaryFromDfCol = memoize(createCategorySummaryFromDfCol);

  getSelectionState(categorySummary) {
    const { categoricalSelection, metadataField } = this.props;
    return Category.getSelectionState(
      categoricalSelection,
      metadataField,
      categorySummary
    );
  }
  
  onSortCategoryLabels = () => {
    const { sortDirection } = this.state;
    if (sortDirection === "descending") {
      this.setState({
        ...this.state,
        sortDirection: "ascending"
      })        
    } else if (sortDirection === "ascending") {
      this.setState({
        ...this.state,
        sortDirection: null
      })
    } else { 
      this.setState({
        ...this.state,
        sortDirection: "descending"
      })
    }
  }
  componentDidUpdate = (prevProps) => {
    const { colors, isExpanded } = this.props;
    const { colorMode, colorAccessor } = colors;
    const { colors: colorsPrev } = prevProps;
    const { colorMode: colorModePrev, colorAccessor: colorAccessorPrev } = colorsPrev;

    const continuousColoring = (
      colorMode === "color by continuous metadata" || 
      colorMode === "color by geneset mean expression" || 
      colorMode =="color by expression"
    )    
    const continuousColoringPrev = (
      colorModePrev === "color by continuous metadata" || 
      colorModePrev === "color by geneset mean expression" || 
      colorModePrev =="color by expression"
    )        
    if (continuousColoring !== continuousColoringPrev && continuousColoringPrev) {
      this.setState({
        ...this.state,
        sortDirection: null
      })
    }

    if ( (colorAccessor !== colorAccessorPrev && isExpanded) || (isExpanded && (isExpanded !== prevProps.isExpanded) && continuousColoring)) {
      const { annoMatrix, metadataField, colors } = this.props;
      this.fetchData(annoMatrix, metadataField, colors).then((val) => {
        const [ categoryData, categorySummary, colorData ] = val;
        if (colorData && categoryData && categorySummary) {
          const col = colorData.col(colorData.colIndex.rindex[0]).asArray();
          const cats = categoryData.col(categoryData.colIndex.rindex[0]).asArray();          

          const averages = {}
          for (const [i,value] of col.entries()){
            averages[cats[i]] = ((averages[cats[i]] ?? 0) + value) / categorySummary.categoryValueCounts[categorySummary.categoryValueIndices.get(cats[i])]
          }
          this.setState({
            continuousAverages: averages
          })
        }
      })
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
  
  handleSankeyClick = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "sankey: toggle",
      category: metadataField,
    });
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

  async fetchData(annoMatrix, metadataField, colors) {
    /*
    fetch our data and the color-by data if appropriate, and then build a summary
    of our category and a color table for the color-by annotation.
    */
    const { schema } = annoMatrix;
    const { colorAccessor, colorMode } = colors;
    const { genesets } = this.props;
    let colorDataPromise = Promise.resolve(null);
    if (colorAccessor) {
      const query = createColorQuery(
        colorMode,
        colorAccessor,
        schema,
        genesets
      );
      if (query) colorDataPromise = annoMatrix.fetch(...query);
    }
    const [categoryData, colorData] = await Promise.all([
      annoMatrix.fetch("obs", metadataField),
      colorDataPromise,
    ]);

    // our data
    const column = categoryData.icol(0);
    const colSchema = schema.annotations.obsByName[metadataField];
    const categorySummary = this.createCategorySummaryFromDfCol(
      column,
      colSchema
    ); 
    return [categoryData, categorySummary, colorData];
  }

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
        categorySummary.allCategoryValues,
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
        categorySummary.allCategoryValues,
        true
      )
    );
  }

  render() {
    const {
      metadataField,
      isExpanded,
      categoricalSelection,
      crossfilter,
      colors,
      annoMatrix,
      isUserAnno,
      sankeySelected,
      layoutChoiceSankey
    } = this.props;
    const { colorMode } = colors;
    const continuousColoring = (colorMode === "color by continuous metadata" || colorMode === "color by geneset mean expression" || colorMode =="color by expression")
    const checkboxID = `category-select-${metadataField}`;
    const sankeyCheckboxID = `sankey-select-${metadataField}`;
    const { sortDirection, continuousAverages } = this.state;

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
                colorAccessor,
                colorTable,
                colorData,
                categoryData,
                categorySummary,
                isColorAccessor,
                handleCategoryToggleAllClick,
              } = asyncProps;
              const isTruncated = !!categorySummary?.isTruncated;
              const selectionState = this.getSelectionState(categorySummary);
              return (
                <CategoryRender
                  metadataField={metadataField}
                  checkboxID={checkboxID}
                  isUserAnno={isUserAnno}
                  isTruncated={isTruncated}
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
                  onCategorySankeyClick={this.handleSankeyClick}
                  sankeyCheckboxID={sankeyCheckboxID}
                  sankeySelected={sankeySelected}
                  layoutChoiceSankey={layoutChoiceSankey}
                  continuousColoring={continuousColoring}
                  sortDirection={sortDirection}
                  onSortCategoryLabels={this.onSortCategoryLabels}
                  continuousAverages={continuousAverages}
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
    sankeyCheckboxID,
    isUserAnno,
    isTruncated,
    isColorAccessor,
    isExpanded,
    selectionState,
    onColorChangeClick,
    onCategoryMenuClick,
    onCategoryMenuKeyPress,
    onCategoryToggleAllClick,
    sankeySelected,
    onCategorySankeyClick,
    layoutChoiceSankey,
    continuousColoring,
    onSortCategoryLabels,
    sortDirection
  }) => {
    /*
    Render category name and controls (eg, color-by button).
    */
    const checkboxRef = useRef(null);
    const checkboxSankeyRef = useRef(null);
    useEffect(() => {
      checkboxRef.current.indeterminate = selectionState === "some";
    }, [checkboxRef.current, selectionState]);
    
    let sortIcon = "expand-all";
    
    if (sortDirection === "ascending"){
      sortIcon="chevron-up"
    } else if (sortDirection === "descending"){
      sortIcon="chevron-down"
    }

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
        {<AnnoDialogEditCategoryName metadataField={metadataField} />}
        {<AnnoDialogAddLabel metadataField={metadataField} />}
        <div>
          <AnchorButton
            onClick={onSortCategoryLabels}
            active={sortDirection}
            minimal
            disabled={!continuousColoring || !isExpanded}
            icon={sortIcon}
          />
          <AnnoMenu
            metadataField={metadataField}
            isUserAnno={isUserAnno}
            createText="Add a new label to this category"
            editText="Edit this category's name"
            deleteText="Delete this category, all associated labels, and remove all cell assignments"
            disableDelete={sankeySelected}
          />

          <Tooltip
            content={
              isTruncated
                ? `Coloring by ${metadataField} is disabled, as it exceeds the limit of ${globals.maxCategoricalOptionsToDisplay} labels`
                : "Use as color scale"
            }
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
              disabled={isTruncated}
              icon="tint"
            />
          </Tooltip>
          <Tooltip
            content="Select category."
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelay}          
          >
            <input
                id={sankeyCheckboxID}
                data-testclass="category-sankey"
                data-testid={`${metadataField}:category-sankey`}
                onChange={onCategorySankeyClick}
                ref={checkboxSankeyRef}
                checked={sankeySelected}
                type="checkbox"
                disabled={layoutChoiceSankey}
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
    sankeySelected,
    onCategorySankeyClick,
    sankeyCheckboxID,
    layoutChoiceSankey,
    continuousColoring,
    onSortCategoryLabels,
    sortDirection,
    continuousAverages
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
            metadataField={metadataField}
            checkboxID={checkboxID}
            sankeyCheckboxID={sankeyCheckboxID}
            isUserAnno={isUserAnno}
            isTruncated={isTruncated}
            isExpanded={isExpanded}
            isColorAccessor={isColorAccessor}
            selectionState={selectionState}
            onColorChangeClick={onColorChangeClick}
            onCategoryToggleAllClick={onCategoryToggleAllClick}
            onCategoryMenuClick={onCategoryMenuClick}
            onCategoryMenuKeyPress={onCategoryMenuKeyPress}
            onCategorySankeyClick={onCategorySankeyClick}
            sankeySelected={sankeySelected}
            layoutChoiceSankey={layoutChoiceSankey}
            continuousColoring={continuousColoring}
            onSortCategoryLabels={onSortCategoryLabels}
            sortDirection={sortDirection}
            isEpxanded={isExpanded}
          />
        </div>
        <div style={{ marginLeft: 26 }}>
          {
            /* values*/
            isExpanded ? (
              <CategoryValueList
                isUserAnno={isUserAnno}
                metadataField={metadataField}
                categoryData={categoryData}
                categorySummary={categorySummary}
                colorAccessor={colorAccessor}
                colorData={colorData}
                colorTable={colorTable}
                sortDirection={sortDirection}
                continuousAverages={continuousAverages}
              />
            ) : null
          }
        </div>
        <div>
          {isExpanded && isTruncated ? (
            <p style={{ paddingLeft: 15 }}>... truncated list ...</p>
          ) : null}
        </div>
      </div>
    );
  }
);

const CategoryValueList = React.memo(
  ({
    isUserAnno,
    metadataField,
    categoryData,
    categorySummary,
    colorAccessor,
    colorData,
    colorTable,
    sortDirection,
    continuousAverages
  }) => {
    let tuples = [...categorySummary.categoryValueIndices];
    
    let newTuples;
    if (continuousAverages) {
      if (sortDirection){
        newTuples = [];
        for (const item in continuousAverages) {
          newTuples.push([item, continuousAverages[item]]);
        }
        
        newTuples.sort(function(a, b) {
            return (a[1] - b[1]);
        });
        for (let i = 0; i < newTuples.length; i++) {
          newTuples[i][1] = i
        }        
        if (sortDirection === "descending"){
          newTuples.reverse()
          for (let i = 0; i < newTuples.length; i++) {
            newTuples[i][1] = i
          }                 
        }
        if (!("unassigned" in continuousAverages)){
          newTuples.push(["unassigned",newTuples.length])
        }
      } else {
        newTuples = tuples;
      }
    } else {
      newTuples = tuples;
    }    
    const newCategoryValues = [];
    const categoryValueCountsObj = {}
    categorySummary.categoryValueCounts.forEach((item,index)=>{
      categoryValueCountsObj[categorySummary.categoryValues[index]] = item;
    })
    const newCategoryValueCounts = []

    for (let i = 0; i < categorySummary.categoryValues.length; i++) {
      newCategoryValues.push(newTuples[i][0])
      newCategoryValueCounts.push(categoryValueCountsObj[newTuples[i][0]])
    }
    const newCategoryValueIndices = new Map(newTuples)
    const newCategorySummary = {
      ...categorySummary, 
      categoryValues: newCategoryValues, 
      categoryValueCounts: newCategoryValueCounts, 
      categoryValueIndices: newCategoryValueIndices
    }    

    if (!isUserAnno) {
      return (
        <>
          {newTuples.map(([value, index]) => (
            <Value
              key={value}
              isUserAnno={isUserAnno}
              metadataField={metadataField}
              categoryIndex={index}
              categoryData={categoryData}
              categorySummary={newCategorySummary}
              colorAccessor={colorAccessor}
              colorData={colorData}
              colorTable={colorTable}
            />
          ))}
        </>
      );
    }

    /* User annotation */
    const flipKey = [...categorySummary.categoryValueIndices].map((t) => t[0]).join("");

    return (
      <Flipper flipKey={flipKey}>
        {newTuples.map(([value, index]) => (
          <Flipped key={value} flipId={value}>
            <Value
              isUserAnno={isUserAnno}
              metadataField={metadataField}
              categoryIndex={index}
              categoryData={categoryData}
              categorySummary={newCategorySummary}
              colorAccessor={colorAccessor}
              colorData={colorData}
              colorTable={colorTable}
              sortDirection={sortDirection}
            />
          </Flipped>
        ))}
      </Flipper>
    );
  }
);
