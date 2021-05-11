/* core dependencies */
import {
  Breadcrumb,
  Menu,
  MenuItem,
  Popover,
  Position,
} from "@blueprintjs/core";
import React from "react";

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
        <Menu style={{ maxHeight: 305, overflow: "scroll" }}>
          {buildDatasetMenuItems(datasets)}
        </Menu>
      }
      position={Position.BOTTOM_LEFT}
    >
      <Breadcrumb>{selectedDatasetName}</Breadcrumb>
    </Popover>
  );
});

export default DatasetMenu;
