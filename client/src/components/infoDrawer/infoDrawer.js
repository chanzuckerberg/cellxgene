import React, { PureComponent } from "react";
import { connect, shallowEqual } from "react-redux";
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
    singleContinuousValues: state.singleContinuousValue.singleContinuousValues,
  };
})
class InfoDrawer extends PureComponent {
  static watchAsync(props, prevProps) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

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
      singleContinuousValues,
    } = this.props;

    const allCategoryNames = selectableCategoryNames(schema).sort();
    const allSingleValues = new Map();

    allCategoryNames.forEach((catName) => {
      const isUserAnno = schema?.annotations?.obsByName[catName]?.writable;
      const colSchema = schema.annotations.obsByName[catName];
      if (!isUserAnno && colSchema.categories?.length === 1) {
        allSingleValues.set(catName, colSchema.categories[0]);
      }
    });
    singleContinuousValues.forEach((value, catName) => {
      allSingleValues.set(catName, value);
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
            allSingleValues,
            dataPortalProps,
          }}
        />
      </Drawer>
    );
  }
}
export default InfoDrawer;
