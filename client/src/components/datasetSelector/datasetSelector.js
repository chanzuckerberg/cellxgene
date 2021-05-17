/* core dependencies */
import { Breadcrumbs } from "@blueprintjs/core";
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

  renderBreadcrumbMenu = ({ datasets, selectedDatasetId }) => {
    /*
    clicking on dataset name opens menu containing all dataset names for the current collection
     */
    const selectedDataset = datasets.find(
      (dataset) => dataset.id === selectedDatasetId
    );
    const datasetsExceptSelected = datasets.filter(
      (dataset) => dataset !== selectedDataset
    );
    return (
      <DatasetMenu
        datasets={datasetsExceptSelected}
        selectedDatasetName={selectedDataset.name}
      />
    );
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
          currentBreadcrumbRenderer={this.renderBreadcrumbMenu}
        />
      </div>
    );
  }
}

export default DatasetSelector;
