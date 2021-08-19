import React, { PureComponent } from "react";
import { connect } from "react-redux";
import { Drawer } from "@blueprintjs/core";

import InfoFormat from "./infoFormat";
import { selectableCategoryNames } from "../../util/stateManager/controlsHelpers";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => {
  return {
    schema: (state as any).annoMatrix.schema,
    datasetTitle: (state as any).config?.displayNames?.dataset ?? "",
    aboutURL: (state as any).config?.links?.["about-dataset"],
    isOpen: (state as any).controls.datasetDrawer,
    dataPortalProps: (state as any).config?.corpora_props,
  };
})
class InfoDrawer extends PureComponent {
  handleClose = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;

    dispatch({ type: "toggle dataset drawer" });
  };

  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'position' does not exist on type 'Readon... Remove this comment to see the full error message
      position,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'aboutURL' does not exist on type 'Readon... Remove this comment to see the full error message
      aboutURL,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'datasetTitle' does not exist on type 'Re... Remove this comment to see the full error message
      datasetTitle,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
      schema,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isOpen' does not exist on type 'Readonly... Remove this comment to see the full error message
      isOpen,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'dataPortalProps' does not exist on type ... Remove this comment to see the full error message
      dataPortalProps,
    } = this.props;

    // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
    const allCategoryNames = selectableCategoryNames(schema).sort();
    const singleValueCategories = new Map();

    allCategoryNames.forEach((catName: any) => {
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
