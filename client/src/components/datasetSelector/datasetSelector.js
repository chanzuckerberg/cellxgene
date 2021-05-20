/* core dependencies */
import { Breadcrumb, Breadcrumbs, Icon } from "@blueprintjs/core";
import React, { PureComponent } from "react";
import { connect } from "react-redux";

/* app dependencies */
import DatasetMenu from "./datasetMenu";

/* styles */
import styles from "./datasetSelector.css";

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
    return { ...breadcrumbProp, className: styles.datasetBreadcrumb };
  };

  buildBreadcrumbProps = (collection, selectedDatasetId) => {
    /*
    create the set of breadcrumbs elements, home > collection name > dataset name, where dataset name reveals the
    dataset menu
     */
    const origin = "https://cellxgene.cziscience.com/"; // TODO(cc) update to ux-dev URL
    const homeProp = this.buildBreadcrumbProp({
      href: origin,
      text: "Home",
    });
    const collectionProp = this.buildBreadcrumbProp({
      href: `${origin}collections/${collection.id}`,
      text: collection.name,
    });
    const datasetProp = this.buildBreadcrumbProp({
      datasets: collection.datasets,
      selectedDatasetId,
    });
    return [homeProp, collectionProp, datasetProp];
  };

  renderBreadcrumb = (selectedDatasetName, renderAsMenu) => {
    /*
    collection name rendered as the last breadcrumb
     */
    const className = renderAsMenu
      ? styles.datasetBreadcrumb
      : styles.datasetDisabledBreadcrumb; /* no sibling datasets */
    return (
      <Breadcrumb className={className}>
        {selectedDatasetName}
        <Icon
          icon="chevron-down"
          style={{ marginLeft: "5px", marginRight: 0 }}
        />
      </Breadcrumb>
    );
  };

  renderBreadcrumbMenu = (datasetsExceptSelected, selectedDatasetName) => {
    /*
    clicking on dataset name opens menu containing all dataset names except the current dataset name for the current
    collection
     */
    return (
      <DatasetMenu datasets={datasetsExceptSelected}>
        {this.renderBreadcrumb(selectedDatasetName, true)}
      </DatasetMenu>
    );
  };

  renderDatasetBreadcrumb = ({ datasets, selectedDatasetId }) => {
    /*
    renders the final dataset breadcrumb; where sibling datasets are selectable by a breadcrumb menu
     */
    const selectedDataset = datasets.find(
      (dataset) => dataset.id === selectedDatasetId
    );
    const siblingDatasets = datasets.filter(
      (dataset) => dataset !== selectedDataset
    );
    const selectedDatasetName = selectedDataset?.name;
    const renderMenu = siblingDatasets.length > 0;
    return renderMenu
      ? this.renderBreadcrumbMenu(siblingDatasets, selectedDatasetName)
      : this.renderBreadcrumb(selectedDatasetName);
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
        }}
      >
        <Breadcrumbs
          items={this.buildBreadcrumbProps(collection, selectedDatasetId)}
          currentBreadcrumbRenderer={this.renderDatasetBreadcrumb}
        />
      </div>
    );
  }
}

export default DatasetSelector;
