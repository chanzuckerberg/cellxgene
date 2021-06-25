/* core dependencies */
import { Breadcrumb, Icon } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { PureComponent } from "react";
import { connect } from "react-redux";

/* app dependencies */
import { switchDataset } from "../../actions";
import DatasetMenu from "./datasetMenu";
import * as globals from "../../globals";
import TruncatingBreadcrumbs from "./truncatingBreadcrumbs";
import { sortDatasets } from "../../util/stateManager/collectionsHelpers";

/* styles */
import styles from "./datasetSelector.css";

/*
app-level collection and dataset breadcrumbs.
 */
@connect((state) => {
  return {
    collection: state.collections?.collection,
    selectedDatasetId: state.collections?.selectedDatasetId,
  };
})
class DatasetSelector extends PureComponent {
  buildBreadcrumbProp = (breadcrumbProp) => {
    /*
    Return base breadcrumb object.
     */
    return { ...breadcrumbProp, className: styles.datasetBreadcrumb };
  };

  buildBreadcrumbProps = (dispatch, collection, selectedDatasetId) => {
    /*
    Create the set of breadcrumbs elements, home > collection name > dataset name, where dataset name reveals the
    dataset menu.
     */
    const origin = globals.API.portalUrl;
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
    const datasets = [...collection.datasets].sort(sortDatasets);
    const datasetProp = this.buildBreadcrumbProp({
      shortText: "Dataset",
      text: selectedDataset.name,
      datasets: datasets.map((dataset) => ({
        ...dataset,
        onClick: () => {
          dispatch(switchDataset(dataset));
        },
      })),
      selectedDatasetId,
    });
    return [homeProp, collectionProp, datasetProp];
  };

  findDatasetById = (selectedDatasetId, datasets) => {
    /*
    Returns the dataset with the given ID.
     */
    return datasets.find((dataset) => dataset.id === selectedDatasetId);
  };

  listSiblingDatasets = (datasets, selectedDataset) => {
    /*
    Returns the set of datasets excluding the given selected dataset.
     */
    return datasets.filter((dataset) => dataset !== selectedDataset);
  };

  renderBreadcrumb = (item, disabled, renderAsMenu) => {
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

  renderBreadcrumbMenu = (item, datasetsExceptSelected) => {
    /*
    Clicking on dataset name opens menu containing all dataset names except the current dataset name for the current
    collection.
     */
    return (
      <DatasetMenu datasets={datasetsExceptSelected}>
        {this.renderBreadcrumb(item, false, true)}
      </DatasetMenu>
    );
  };

  renderBreadcrumbMenuIcon = (renderAsMenu) => {
    /*
    Render breadcrumb menu icon "chevron down".
     */
    return renderAsMenu ? (
      <Icon
        icon={IconNames.CHEVRON_DOWN}
        style={{ marginLeft: "5px", marginRight: 0 }}
      />
    ) : null;
  };

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

  render() {
    const { collection, dispatch, selectedDatasetId } = this.props;
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
          breadcrumbRenderer={this.renderBreadcrumb}
          currentBreadcrumbRenderer={this.renderDatasetBreadcrumb}
          items={this.buildBreadcrumbProps(
            dispatch,
            collection,
            selectedDatasetId
          )}
        />
      </div>
    );
  }
}

export default DatasetSelector;
