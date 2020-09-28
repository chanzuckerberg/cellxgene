import React, { PureComponent } from "react";
import { connect, shallowEqual } from "react-redux";
import { Drawer } from "@blueprintjs/core";
import Async from "react-async";

import InfoFormat from "./infoFormat";
import {
  selectableCategoryNames,
  createCategorySummaryFromDfCol,
} from "../../util/stateManager/controlsHelpers";

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

    return (
      <Drawer
        title="Dataset Overview"
        onClose={this.handleClose}
        {...{ isOpen, position }}
      >
        <Async
          watchFn={InfoDrawer.watchAsync}
          promiseFn={this.fetchAsyncProps}
          watchProps={{ schema }}
        >
          <Async.Pending>
            <InfoFormat
              skeleton
              {...{ datasetTitle, aboutURL, dataPortalProps }}
            />
          </Async.Pending>
          <Async.Rejected>
            {(error) => {
              console.error(error);
              return <span>Failed to load info</span>;
            }}
          </Async.Rejected>
          <Async.Fulfilled>
            {(asyncProps) => {
              const { singleValueCategories } = asyncProps;
              return (
                <InfoFormat
                  {...{
                    datasetTitle,
                    aboutURL,
                    singleValueCategories,
                    dataPortalProps,
                  }}
                />
              );
            }}
          </Async.Fulfilled>
        </Async>
      </Drawer>
    );
  }
}
export default InfoDrawer;
