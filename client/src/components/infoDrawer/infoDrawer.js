import React, { PureComponent } from "react";
import { connect, shallowEqual } from "react-redux";
import { Drawer } from "@blueprintjs/core";
import Async from "react-async";
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

    const obsIndex = schema.annotations.obs.index;
    const allContinuousNames = schema.annotations.obs.columns
      .filter((col) => col.type === "int32" || col.type === "float32")
      .filter((col) => col.name !== obsIndex)
      .filter((col) => !col.writable) // skip user annotations - they will be treated as categorical
      .map((col) => col.name);
    const annoContinous = allContinuousNames.map((catName) => {
      return annoMatrix.fetch("obs", catName);
    });
    const singleValueContinuous = (await Promise.all(annoContinous)).reduce(
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
    return { singleValueContinuous };
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
        <Async
          watchFn={InfoDrawer.watchAsync}
          promiseFn={this.fetchAsyncProps}
          watchProps={{ schema }}
        >
          <Async.Pending>
            <InfoFormat
              {...{
                datasetTitle,
                aboutURL,
                singleValueCategories,
                dataPortalProps,
              }}
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
              const { singleValueContinuous } = asyncProps;
              for (const [key, value] of singleValueContinuous.entries()) {
                singleValueCategories.set(key, value);
              }
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
