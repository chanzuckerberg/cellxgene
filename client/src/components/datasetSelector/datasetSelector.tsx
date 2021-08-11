/* core dependencies */
import { Breadcrumb, Icon } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { PureComponent } from "react";
import { connect } from "react-redux";

/* app dependencies */
import { openDataset, switchDataset } from "../../actions";
import DatasetMenu from "./datasetMenu";
import * as globals from "../../globals";
import TruncatingBreadcrumbs from "./truncatingBreadcrumbs";
import { sortDatasets } from "../../util/stateManager/collectionsHelpers";

/* styles */
// @ts-expect-error --- TODO revisit
import styles from "./datasetSelector.css";

/*
app-level collection and dataset breadcrumbs.
 */
// @ts-expect-error ts-migrate(1238) TODO revisit
@connect((state) => {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- TODO revisit
  const genesetsInProgress = (state as any).genesets?.genesets?.size > 0;
  const individualGenesInProgress =
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- TODO revisit
    (state as any).controls?.userDefinedGenes?.length > 0;
  return {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- TODO revisit
    collection: (state as any).collections?.collection,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- TODO revisit
    selectedDatasetId: (state as any).collections?.selectedDatasetId,
    workInProgress: genesetsInProgress || individualGenesInProgress,
  };
})
class DatasetSelector extends PureComponent {
  // @ts-expect-error --- TODO revisit
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
  buildBreadcrumbProp = (breadcrumbProp) =>
    /*
    Return base breadcrumb object.
     */
    ({ ...breadcrumbProp, className: styles.datasetBreadcrumb });

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
  buildBreadcrumbProps = (
    // @ts-expect-error --- TODO revisit
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
    dispatch,
    // @ts-expect-error --- TODO revisit
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
    collection,
    // @ts-expect-error --- TODO revisit
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
    selectedDatasetId,
    // @ts-expect-error --- TODO revisit
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
    workInProgress
  ) => {
    /*
    Create the set of breadcrumbs elements, home > collection name > dataset name, where dataset name reveals the
    dataset menu.
     */
    const { origin } = globals.API;
    const homeProp = this.buildBreadcrumbProp({
      href: origin,
      shortText: "Home",
      text: "Home",
    });
    const collectionProp = this.buildBreadcrumbProp({
      href: `${origin}collections/${collection.id}`,
      shortText: "Collection",
      text: collection.name,
    });
    const selectedDataset = this.findDatasetById(
      selectedDatasetId,
      collection.datasets
    );
    const datasets = [...collection.datasets]
      .sort(sortDatasets)
      .map((dataset) => {
        const dispatchAction = workInProgress
          ? openDataset(dataset)
          : switchDataset(dataset);
        return {
          ...dataset,
          onClick: () => {
            dispatch(dispatchAction);
          },
        };
      });
    const datasetProp = this.buildBreadcrumbProp({
      shortText: "Dataset",
      text: selectedDataset.name,
      datasets,
      selectedDatasetId,
    });
    return [homeProp, collectionProp, datasetProp];
  };

  /*
  Returns the dataset with the given ID.
   */
  // @ts-expect-error --- TODO revisit
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
  findDatasetById = (selectedDatasetId, datasets) =>
    // @ts-expect-error --- TODO revisit
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
    datasets.find((dataset) => dataset.id === selectedDatasetId);

  /*
   Returns the set of datasets excluding the given selected dataset.
   */
  // @ts-expect-error --- TODO revisit
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
  listSiblingDatasets = (datasets, selectedDataset) =>
    // @ts-expect-error --- TODO revisit
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
    datasets.filter((dataset) => dataset !== selectedDataset);

  // @ts-expect-error --- TODO revisit
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
  renderBreadcrumb = (item, disabled, renderAsMenu?) => {
    /*
    Render BP Breadcrumb, adding menu-specific styles if necessary.
    TODO(cc) split and simplify breadcrumb versus menu breadcrumb functionality.
     */
    const className = disabled
      ? styles.datasetDisabledBreadcrumb /* no sibling datasets */
      : styles.datasetBreadcrumb;
    return (
      <Breadcrumb href={item.href} className={className}>
        {item.displayText}
        {this.renderBreadcrumbMenuIcon(renderAsMenu)}
      </Breadcrumb>
    );
  };

  // @ts-expect-error --- TODO revisit
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
  renderBreadcrumbMenu = (item, datasetsExceptSelected) => 
    /*
    Clicking on dataset name opens menu containing all dataset names except the current dataset name for the current
    collection.
     */
     (
      // @ts-expect-error --- TODO revisit
      <DatasetMenu datasets={datasetsExceptSelected}>
        {this.renderBreadcrumb(item, false, true)}
      </DatasetMenu>
    )
  ;

  /*
   Render breadcrumb menu icon "chevron down".
   */
  // @ts-expect-error --- TODO revisit
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
  renderBreadcrumbMenuIcon = (renderAsMenu) =>
    renderAsMenu ? (
      <Icon
        icon={IconNames.CHEVRON_DOWN}
        style={{ marginLeft: "5px", marginRight: 0 }}
      />
    ) : null;

  // @ts-expect-error --- TODO revisit
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
  renderDatasetBreadcrumb = (item) => {
    /*
     Renders the final dataset breadcrumb where sibling datasets are selectable by a breadcrumb menu.
     */
    const { datasets, selectedDatasetId } = item;
    const selectedDataset = this.findDatasetById(selectedDatasetId, datasets);
    const siblingDatasets = this.listSiblingDatasets(datasets, selectedDataset);
    const renderMenu = siblingDatasets.length > 0;
    return renderMenu
      ? this.renderBreadcrumbMenu(item, siblingDatasets)
      : this.renderBreadcrumb(item, true);
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO add return value
  render() {
    // @ts-expect-error --- TODO revisit
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
    const { collection, dispatch, selectedDatasetId, workInProgress } =
      this.props;
    if (!collection) {
      return null;
    }
    return (
      <div
        style={{
          marginTop: "8px", // Match margin on sibling menu buttons
          flexGrow: 1,
          overflow: "scroll",
          flex: 1,
        }}
      >
        <TruncatingBreadcrumbs
          // @ts-expect-error --- TODO revisit
          // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- TODO revisit
          breadcrumbRenderer={this.renderBreadcrumb}
          currentBreadcrumbRenderer={this.renderDatasetBreadcrumb}
          items={this.buildBreadcrumbProps(
            dispatch,
            collection,
            selectedDatasetId,
            workInProgress
          )}
        />
      </div>
    );
  }
}

export default DatasetSelector;
