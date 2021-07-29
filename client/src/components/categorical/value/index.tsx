import { connect } from "react-redux";
import React from "react";
import * as d3 from "d3";

import {
  Button,
  Classes,
  Icon,
  Menu,
  MenuItem,
  Popover,
  PopoverInteractionKind,
  Position,
} from "@blueprintjs/core";
import * as globals from "../../../globals";
// @ts-expect-error ts-migrate(2307) FIXME: Cannot find module '../categorical.css' or its cor... Remove this comment to see the full error message
import styles from "../categorical.css";
import AnnoDialog from "../../annoDialog";
import LabelInput from "../../labelInput";
import Truncate from "../../util/truncate";

import { AnnotationsHelpers } from "../../../util/stateManager";
import { labelPrompt, isLabelErroneous } from "../labelUtil";
import actions from "../../../actions";
import MiniHistogram from "../../miniHistogram";
import MiniStackedBar from "../../miniStackedBar";
import { CategoryCrossfilterContext } from "../categoryContext";

const STACKED_BAR_HEIGHT = 11;
const STACKED_BAR_WIDTH = 100;

/* this is defined outside of the class so we can use it in connect() */
// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function _currentLabelAsString(ownProps: any) {
  const { label } = ownProps;
  // when called as a function, the String() constructor performs type conversion,
  // and returns a primitive string.
  return String(label);
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state, ownProps) => {
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'pointDilation' does not exist on type 'D... Remove this comment to see the full error message
  const { pointDilation, categoricalSelection } = state;
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type '{... Remove this comment to see the full error message
  const { metadataField, categorySummary, categoryIndex } = ownProps;
  const isDilated =
    pointDilation.metadataField === metadataField &&
    pointDilation.categoryField === _currentLabelAsString(ownProps);

  const category = categoricalSelection[metadataField];
  const label = categorySummary.categoryValues[categoryIndex];
  const isSelected = category.get(label) ?? true;

  return {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    annotations: (state as any).annotations,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    schema: (state as any).annoMatrix?.schema,
    isDilated,
    isSelected,
    label,
  };
})
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class CategoryValue extends React.Component<{}, State> {
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      editedLabelText: this.currentLabelAsString(),
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  componentDidUpdate(prevProps: {}) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
    const { metadataField, categoryIndex, categorySummary } = this.props;
    if (
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (prevProps as any).metadataField !== metadataField ||
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (prevProps as any).categoryIndex !== categoryIndex || // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (prevProps as any).categorySummary !== categorySummary
    ) {
      // eslint-disable-next-line react/no-did-update-set-state --- adequately checked to prevent looping
      this.setState({
        editedLabelText: this.currentLabelAsString(),
      });
    }
  }

  // If coloring by and this isn't the colorAccessor and it isn't being edited
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  get shouldRenderStackedBarOrHistogram() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorAccessor' does not exist on type 'R... Remove this comment to see the full error message
    const { colorAccessor, isColorBy, annotations } = this.props;

    return !!colorAccessor && !isColorBy && !annotations.isEditingLabelName;
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleDeleteValue = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField, label } = this.props;
    dispatch(actions.annotationDeleteLabelFromCategory(metadataField, label));
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleAddCurrentSelectionToThisLabel = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField, label } = this.props;
    dispatch(actions.annotationLabelCurrentSelection(metadataField, label));
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleEditValue = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField, label } = this.props;
    const { editedLabelText } = this.state;
    this.cancelEditMode();
    dispatch(
      actions.annotationRenameLabelInCategory(
        metadataField,
        label,
        editedLabelText
      )
    );
    e.preventDefault();
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleCreateArbitraryLabel = (txt: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField, label } = this.props;
    this.cancelEditMode();
    dispatch(
      actions.annotationRenameLabelInCategory(metadataField, label, txt)
    );
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  labelNameError = (name: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
    const { metadataField, schema } = this.props;
    if (name === this.currentLabelAsString()) return false;
    return isLabelErroneous(name, metadataField, schema);
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  instruction = (label: any) => {
    return labelPrompt(this.labelNameError(label), "New, unique label", ":");
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  activateEditLabelMode = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField, categoryIndex, label } = this.props;
    dispatch({
      type: "annotation: activate edit label mode",
      metadataField,
      categoryIndex,
      label,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  cancelEditMode = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField, categoryIndex, label } = this.props;
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  toggleOff = () => {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
      dispatch,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
      metadataField,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryIndex' does not exist on type 'R... Remove this comment to see the full error message
      categoryIndex,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categorySummary' does not exist on type ... Remove this comment to see the full error message
      categorySummary,
    } = this.props;
    const label = categorySummary.categoryValues[categoryIndex];
    dispatch(
      actions.selectCategoricalMetadataAction(
        "categorical metadata filter deselect",
        metadataField,
        categorySummary.allCategoryValues,
        label,
        false
      )
    );
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  shouldComponentUpdate = (nextProps: any, nextState: any) => {
    /*
    Checks to see if at least one of the following changed:
    * world state
    * the color accessor (what is currently being colored by)
    * if this categorical value's selection status has changed
    * the crossfilter (ie, global selection state)

    If and only if true, update the component
    */
    const { props, state } = this;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryIndex' does not exist on type 'R... Remove this comment to see the full error message
    const { categoryIndex, categorySummary, isSelected } = props;
    const {
      categoryIndex: newCategoryIndex,
      categorySummary: newCategorySummary,
      isSelected: newIsSelected,
    } = nextProps;

    const label = categorySummary.categoryValues[categoryIndex];
    const newLabel = newCategorySummary.categoryValues[newCategoryIndex];
    const labelChanged = label !== newLabel;
    const valueSelectionChange = isSelected !== newIsSelected;

    const colorAccessorChange =
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (props as any).colorAccessor !== nextProps.colorAccessor;
    const annotationsChange =
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (props as any).annotations !== nextProps.annotations;
    const editingLabel = state.editedLabelText !== nextState.editedLabelText;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const dilationChange = (props as any).isDilated !== nextProps.isDilated;

    const count = categorySummary.categoryValueCounts[categoryIndex];
    const newCount = newCategorySummary.categoryValueCounts[newCategoryIndex];
    const countChanged = count !== newCount;

    // If the user edits an annotation that is currently colored-by, colors may be re-assigned.
    // This test is conservative - it may cause re-rendering of entire category (all labels)
    // if any one changes, but only for the currently colored-by category.
    const colorMightHaveChanged =
      nextProps.colorAccessor === nextProps.metadataField &&
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (props as any).categorySummary !== nextProps.categorySummary;

    return (
      labelChanged ||
      valueSelectionChange ||
      colorAccessorChange ||
      annotationsChange ||
      editingLabel ||
      dilationChange ||
      countChanged ||
      colorMightHaveChanged
    );
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  toggleOn = () => {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
      dispatch,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
      metadataField,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryIndex' does not exist on type 'R... Remove this comment to see the full error message
      categoryIndex,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categorySummary' does not exist on type ... Remove this comment to see the full error message
      categorySummary,
    } = this.props;
    const label = categorySummary.categoryValues[categoryIndex];
    dispatch(
      actions.selectCategoricalMetadataAction(
        "categorical metadata filter select",
        metadataField,
        categorySummary.allCategoryValues,
        label,
        true
      )
    );
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleMouseEnter = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField, categoryIndex, label } = this.props;
    dispatch({
      type: "category value mouse hover start",
      metadataField,
      categoryIndex,
      label,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleMouseExit = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField, categoryIndex, label } = this.props;
    dispatch({
      type: "category value mouse hover end",
      metadataField,
      categoryIndex,
      label,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleTextChange = (text: any) => {
    this.setState({ editedLabelText: text });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleChoice = (e: any) => {
    /* Blueprint Suggest format */
    this.setState({ editedLabelText: e.target });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  createHistogramBins = (
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    metadataField: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    categoryData: any,
    // @ts-expect-error ts-migrate(6133) FIXME: 'colorAccessor' is declared but its value is never... Remove this comment to see the full error message
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    colorAccessor: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    colorData: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    categoryValue: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    width: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    height: any
  ) => {
    /*
      Knowing that colorScale is based off continuous data,
      createHistogramBins fetches the continuous data in relation to the cells relevant to the category value.
      It then separates that data into 50 bins for drawing the mini-histogram
    */
    const groupBy = categoryData.col(metadataField);
    const col = colorData.icol(0);
    const range = col.summarize();

    const histogramMap = col.histogram(
      50,
      [range.min, range.max],
      groupBy
    ); /* Because the signature changes we really need different names for histogram to differentiate signatures  */

    const bins = histogramMap.has(categoryValue)
      ? histogramMap.get(categoryValue)
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  createStackedGraphBins = (
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    metadataField: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    categoryData: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    colorAccessor: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    colorData: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    categoryValue: any,
    // @ts-expect-error ts-migrate(6133) FIXME: 'colorTable' is declared but its value is never re... Remove this comment to see the full error message
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    colorTable: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    schema: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    width: any
  ) => {
    /*
      Knowing that the color scale is based off of categorical data,
      createOccupancyStack obtains a map showing the number if cells per colored value
      Using the colorScale a stack of colored bars is drawn representing the map
     */
    const groupBy = categoryData.col(metadataField);
    const occupancyMap = colorData
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

      const dfColumn = colorData.col(colorAccessor);
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  currentLabelAsString() {
    return _currentLabelAsString(this.props);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  isAddCurrentSelectionDisabled(crossfilter: any, category: any, value: any) {
    /*
    disable "add current selection to label", if one of the following is true:
    1. no cells are selected
    2. all currently selected cells already have this label, on this category
    */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryData' does not exist on type 'Re... Remove this comment to see the full error message
    const { categoryData } = this.props;

    // 1. no cells selected?
    if (crossfilter.countSelected() === 0) {
      return true;
    }
    // 2. all selected cells already have the label
    const mask = crossfilter.allSelectedMask();
    if (
      AnnotationsHelpers.allHaveLabelByMask(categoryData, category, value, mask)
    ) {
      return true;
    }
    // else, don't disable
    return false;
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  renderMiniStackedBar = () => {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorAccessor' does not exist on type 'R... Remove this comment to see the full error message
      colorAccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
      metadataField,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryData' does not exist on type 'Re... Remove this comment to see the full error message
      categoryData,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorData' does not exist on type 'Reado... Remove this comment to see the full error message
      colorData,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorTable' does not exist on type 'Read... Remove this comment to see the full error message
      colorTable,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
      schema,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'label' does not exist on type 'Readonly<... Remove this comment to see the full error message
      label,
    } = this.props;
    const isColorBy = metadataField === colorAccessor;

    if (
      !this.shouldRenderStackedBarOrHistogram ||
      !AnnotationsHelpers.isCategoricalAnnotation(schema, colorAccessor) ||
      isColorBy
    ) {
      return null;
    }

    const { domainValues, scale, domain, occupancy } =
      this.createStackedGraphBins(
        metadataField,
        categoryData,
        colorAccessor,
        colorData,
        label,
        colorTable,
        schema,
        STACKED_BAR_WIDTH
      ) ?? {};

    if (!domainValues || !scale || !domain || !occupancy) {
      return null;
    }

    return (
      <MiniStackedBar
        {...{
          colorTable,
          domainValues,
          scale,
          domain,
          occupancy,
        }}
        // @ts-expect-error ts-migrate(2322) FIXME: Type '{ height: number; width: number; colorTable:... Remove this comment to see the full error message
        height={STACKED_BAR_HEIGHT}
        width={STACKED_BAR_WIDTH}
      />
    );
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  renderMiniHistogram = () => {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorAccessor' does not exist on type 'R... Remove this comment to see the full error message
      colorAccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
      metadataField,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorData' does not exist on type 'Reado... Remove this comment to see the full error message
      colorData,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryData' does not exist on type 'Re... Remove this comment to see the full error message
      categoryData,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorTable' does not exist on type 'Read... Remove this comment to see the full error message
      colorTable,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
      schema,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'label' does not exist on type 'Readonly<... Remove this comment to see the full error message
      label,
    } = this.props;
    const colorScale = colorTable?.scale;

    if (
      !this.shouldRenderStackedBarOrHistogram ||
      // This function returns true on categorical annotations(when stacked bar should not render),
      //  in cases where the colorAccessor is a gene this function will return undefined since genes do not live on the schema
      AnnotationsHelpers.isCategoricalAnnotation(schema, colorAccessor) === true
    ) {
      return null;
    }

    const { xScale, yScale, bins } =
      this.createHistogramBins(
        metadataField,
        categoryData,
        colorAccessor,
        colorData,
        label,
        STACKED_BAR_WIDTH,
        STACKED_BAR_HEIGHT
      ) ?? {}; // if createHistogramBins returns empty object assign null to deconstructed

    if (!xScale || !yScale || !bins) return null;

    return (
      <MiniHistogram
        {...{
          colorScale,
          xScale,
          yScale,
          bins,
        }}
        // @ts-expect-error ts-migrate(2322) FIXME: Type '{ obsOrVarContinuousFieldDisplayName: any; d... Remove this comment to see the full error message
        obsOrVarContinuousFieldDisplayName={colorAccessor}
        domainLabel={label}
        height={STACKED_BAR_HEIGHT}
        width={STACKED_BAR_WIDTH}
      />
    );
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
      metadataField,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryIndex' does not exist on type 'R... Remove this comment to see the full error message
      categoryIndex,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorAccessor' does not exist on type 'R... Remove this comment to see the full error message
      colorAccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorTable' does not exist on type 'Read... Remove this comment to see the full error message
      colorTable,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isUserAnno' does not exist on type 'Read... Remove this comment to see the full error message
      isUserAnno,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'annotations' does not exist on type 'Rea... Remove this comment to see the full error message
      annotations,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isDilated' does not exist on type 'Reado... Remove this comment to see the full error message
      isDilated,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isSelected' does not exist on type 'Read... Remove this comment to see the full error message
      isSelected,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categorySummary' does not exist on type ... Remove this comment to see the full error message
      categorySummary,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'label' does not exist on type 'Readonly<... Remove this comment to see the full error message
      label,
    } = this.props;
    const colorScale = colorTable?.scale;

    const { editedLabelText } = this.state;

    const count = categorySummary.categoryValueCounts[categoryIndex];
    const displayString = this.currentLabelAsString();

    /* this is the color scale, so add swatches below */
    const isColorBy = metadataField === colorAccessor;
    const { categoryValueIndices } = categorySummary;

    const editModeActive =
      isUserAnno &&
      annotations.labelEditable.category === metadataField &&
      annotations.isEditingLabelName &&
      annotations.labelEditable.label === categoryIndex;

    const valueToggleLabel = `value-toggle-checkbox-${metadataField}-${displayString}`;

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
          STACKED_BAR_WIDTH -
          CHART_MARGIN
        : globals.leftSidebarWidth - otherElementsWidth;

    return (
      <div
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
              className={`${Classes.CONTROL} ${Classes.CHECKBOX}`}
              style={{ margin: 0 }}
            >
              <input
                id={valueToggleLabel}
                onChange={isSelected ? this.toggleOff : this.toggleOn}
                data-testclass="categorical-value-select"
                data-testid={`categorical-value-select-${metadataField}-${displayString}`}
                checked={isSelected}
                type="checkbox"
              />
              <span
                className={Classes.CONTROL_INDICATOR}
                onMouseEnter={this.handleMouseExit}
                onMouseLeave={this.handleMouseEnter}
              />
            </label>
            <Truncate>
              <span
                data-testid={`categorical-value-${metadataField}-${displayString}`}
                data-testclass="categorical-value"
                // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'number | ... Remove this comment to see the full error message
                tabIndex="-1"
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
                  // @ts-expect-error ts-migrate(2322) FIXME: Type '{ isActive: any; inputProps: { "data-testid"... Remove this comment to see the full error message
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
                      // @ts-expect-error ts-migrate(2322) FIXME: Type '{ label: any; labelSuggestions: null; onChan... Remove this comment to see the full error message
                      label={editedLabelText}
                      labelSuggestions={null}
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
            {this.renderMiniStackedBar()}
            {this.renderMiniHistogram()}
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
              display={isColorBy && categoryValueIndices ? "auto" : "none"}
              style={{
                marginLeft: 5,
                width: 15,
                height: 15,
                backgroundColor:
                  isColorBy && categoryValueIndices
                    ? colorScale(categoryValueIndices.get(label))
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
                      <CategoryCrossfilterContext.Consumer>
                        {(crossfilter) => (
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
                              crossfilter,
                              metadataField,
                              label
                            )}
                          />
                        )}
                      </CategoryCrossfilterContext.Consumer>
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
                          icon="trash"
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
