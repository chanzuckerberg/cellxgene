import React, { PureComponent } from "react";
import { connect } from "react-redux";
import { Drawer } from "@blueprintjs/core";

import InfoFormat from "./infoFormat";
import { selectableCategoryNames } from "../../util/stateManager/controlsHelpers";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  collection: ((state as any).collections as any)?.collection,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  schema: (state as any).annoMatrix.schema,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  datasetTitle: (state as any).config?.displayNames?.dataset ?? "",
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  aboutURL: (state as any).config?.links?.["about-dataset"],
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  isOpen: (state as any).controls.datasetDrawer,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  dataPortalProps: (state as any).config?.corpora_props,
}))
class InfoDrawer extends PureComponent {
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleClose = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;

    dispatch({ type: "toggle dataset drawer" });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'collection' does not exist on type 'Readon... Remove this comment to see the full error message
      collection,
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

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    allCategoryNames.forEach((catName: any) => {
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
