import React, { PureComponent } from "react";
import { connect } from "react-redux";
import { Drawer } from "@blueprintjs/core";

import InfoFormat from "./infoFormat";
import { selectableCategoryNames } from "../../util/stateManager/controlsHelpers";

@connect((state) => {
  return {
    schema: state.annoMatrix.schema,
    datasetTitle: state.config?.displayNames?.dataset ?? "",
    aboutURL: state.config?.links?.["about-dataset"],
    isOpen: state.controls.datasetDrawer,
    dataPortalProps: state.config?.["corpora_props"],
  };
})
class InfoDrawer extends PureComponent {
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
            dataPortalProps: dataPortalProps ?? {},
          }}
        />
      </Drawer>
    );
  }
}
export default InfoDrawer;
