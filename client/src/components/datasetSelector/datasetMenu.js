/* core dependencies */
import {
  Breadcrumb,
  Colors,
  Icon,
  Menu,
  MenuItem,
  Popover,
  Position,
} from "@blueprintjs/core";
import React from "react";

/* styles */
import styles from "./datasetSelector.css";

const buildDatasetMenuItems = (datasets) => {
  /*
  map dataset to menu item
   */
  return datasets.map((dataset) => {
    return <MenuItem key={dataset.id} href={dataset.url} text={dataset.name} />;
  });
};
/*
 dataset menu, toggled from dataset name in app-level breadcrumbs
 */
const DatasetMenu = React.memo(({ datasets, selectedDatasetName }) => {
  return (
    <Popover
      content={
        <Menu style={{ color: Colors.BLACK }}>
          {buildDatasetMenuItems(datasets)}
        </Menu>
      }
      minimal
      modifiers={{ offset: { offset: "0, 10" } }}
      popoverClassName={styles.datasetPopover}
      position={Position.BOTTOM_LEFT}
      targetClassName={styles.datasetPopoverTarget}
    >
      <Breadcrumb className={styles.datasetBreadcrumb}>
        {selectedDatasetName}
        <Icon
          icon="chevron-down"
          style={{ marginLeft: "5px", marginRight: 0 }}
        />
      </Breadcrumb>
    </Popover>
  );
});
export default DatasetMenu;
