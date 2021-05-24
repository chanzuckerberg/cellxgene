import React, { PureComponent } from "react";
import { connect } from "react-redux";
import { Drawer } from "@blueprintjs/core";

import InfoFormat from "./infoFormat";
import { selectableCategoryNames } from "../../util/stateManager/controlsHelpers";

@connect((state) => {
  const selectedDatasetId = state.collections?.selectedDatasetId;
  const collection = state.collections?.collectionsByDatasetId?.get(
    selectedDatasetId
  );
  return {
    collection,
    schema: state.annoMatrix.schema,
    isOpen: state.controls.datasetDrawer,
    dataPortalProps: state.config?.corpora_props,
  };
})
class InfoDrawer extends PureComponent {
  handleClose = () => {
    const { dispatch } = this.props;

    dispatch({ type: "toggle dataset drawer" });
  };

  render() {
    const {
      collection,
      position,
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
      <Drawer onClose={this.handleClose} size={480} {...{ isOpen, position }}>
        <InfoFormat
          {...{
            collection,
            singleValueCategories,
            dataPortalProps: dataPortalProps ?? {},
          }}
        />
      </Drawer>
    );
  }
}
export default InfoDrawer;
