import React, { PureComponent } from "react";
import { connect } from "react-redux";
import { Drawer } from "@blueprintjs/core";

import InfoFormat from "./infoFormat";
import { selectableCategoryNames } from "../../util/stateManager/controlsHelpers";

@connect((state) => {
  return {
    annoMatrix: state.annoMatrix,
    schema: state.annoMatrix.schema,
    datasetTitle: state.config?.displayNames?.dataset ?? "",
    aboutURL: state.config?.links?.["about-dataset"],
    isOpen: state.controls.datasetDrawer,
    dataPortalProps: state.config?.["corpora_props"] ?? {},
  };
})
class InfoDrawer extends PureComponent {
  static watchAsync(props, prevProps) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  fetchAsyncProps = async (props) => {
    const { schema } = props.watchProps;
    const { annoMatrix } = this.props;

    const allCategoryNames = selectableCategoryNames(schema).sort();
    const obsIndex = schema.annotations.obs.index;
    const allContinuousNames = schema.annotations.obs.columns
      .filter((col) => col.type === "int32" || col.type === "float32")
      .filter((col) => col.name !== obsIndex)
      .filter((col) => !col.writable) // skip user annotations - they will be treated as categorical
      .map((col) => col.name);
    const annoContinous = allContinuousNames.map((catName) => {
      return annoMatrix.fetch("obs", catName);
    });
    const singleValueContinous = (await Promise.all(annoContinous)).reduce(
      (acc, continuousData, i) => {
        if (!continuousData) return acc;
        const column = continuousData.icol(0);
        const catName = allContinuousNames[i];
        const summary = column.summarize();
        if (summary.min === summary.max) {
          acc.set(catName, summary.min);
        }
        return acc;
      },
      new Map()
    );
    const nonUserAnnoCategories = allCategoryNames.map((catName) => {
      const isUserAnno = schema?.annotations?.obsByName[catName]?.writable;
      if (!isUserAnno) return annoMatrix.fetch("obs", catName);
      return null;
    });
    const singleValueCategories = (
      await Promise.all(nonUserAnnoCategories)
    ).reduce((acc, categoryData, i) => {
      // Actually check to see if it is null(user anno)
      if (!categoryData) return acc;
      const catName = allCategoryNames[i];

      const column = categoryData.icol(0);
      const colSchema = schema.annotations.obsByName[catName];

      const categorySummary = createCategorySummaryFromDfCol(column, colSchema);

      const { numCategoryValues } = categorySummary;
      //  Add to the array if the category has only one value
      if (numCategoryValues === 1) {
        acc.set(catName, categorySummary.allCategoryValues[0]);
      }
      return acc;
    }, new Map());
    for (const [key, value] of singleValueContinous.entries()) {
      singleValueCategories.set(key, value);
    }
    return { singleValueCategories };
  };

  handleClose = () => {
    const { dispatch } = this.props;

    dispatch({ type: "toggle dataset drawer" });
  };

  render() {
    const {
      position,
      aboutURL,
      datasetTitle,
      schema,
      isOpen,
      dataPortalProps,
    } = this.props;

    const allCategoryNames = selectableCategoryNames(schema).sort();
    const singleValueCategories = new Map();

    allCategoryNames.forEach((catName) => {
      const isUserAnno = schema?.annotations?.obsByName[catName]?.writable;
      const colSchema = schema.annotations.obsByName[catName];
      if (!isUserAnno && colSchema.categories?.length === 1) {
        singleValueCategories.set(catName, colSchema.categories[0]);
      }
    });

    return (
      <Drawer
        title="Dataset Overview"
        onClose={this.handleClose}
        {...{ isOpen, position }}
      >
        <InfoFormat
          {...{
            datasetTitle,
            aboutURL,
            singleValueCategories,
            dataPortalProps,
          }}
        />
      </Drawer>
    );
  }
}
export default InfoDrawer;
