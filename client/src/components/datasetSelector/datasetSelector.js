/* core dependencies */
import { Breadcrumb, Icon } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { PureComponent } from "react";
import { connect } from "react-redux";

/* app dependencies */
import DatasetMenu from "./datasetMenu";

/* styles */
import styles from "./datasetSelector.css";
import TruncatingBreadcrumbs from "./truncatingBreadcrumbs";

/*
app-level collection and dataset breadcrumbs.
 */
@connect((state) => {
  const selectedDatasetId = state.collections?.selectedDatasetId;
  const collection = state.collections?.collectionsByDatasetId?.get(
    selectedDatasetId
  );
  return {
    collection,
    selectedDatasetId,
  };
})
class DatasetSelector extends PureComponent {
  buildBreadcrumbProp = (breadcrumbProp) => {
    /*
    Return base breadcrumb object.
     */
    return { ...breadcrumbProp, className: styles.datasetBreadcrumb };
  };

  buildBreadcrumbProps = (collection, selectedDatasetId) => {
    /*
    Create the set of breadcrumbs elements, home > collection name > dataset name, where dataset name reveals the
    dataset menu.
     */
    const origin = "https://cellxgene.cziscience.com/"; // TODO(cc) update to ux-dev URL
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
    const datasetProp = this.buildBreadcrumbProp({
      shortText: "Dataset",
      text: selectedDataset.name,
      datasets: collection.datasets,
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
    const { collection, selectedDatasetId } = this.props;
    if (!collection) {
      return null;
    }
    return (
      <div
        style={{
          marginTop: "8px", // Match margin on sibling menu buttons
          flexGrow: 1,
          overflow: "scroll", // TODO(cc) Mim to revisit
          flex: 1,
        }}
      >
        <TruncatingBreadcrumbs
          breadcrumbRenderer={this.renderBreadcrumb}
          currentBreadcrumbRenderer={this.renderDatasetBreadcrumb}
          items={this.buildBreadcrumbProps(collection, selectedDatasetId)}
        />
      </div>
    );
  }
}

export default DatasetSelector;
