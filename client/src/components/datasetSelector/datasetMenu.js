/* core dependencies */
import { Menu, MenuItem, Popover, Position } from "@blueprintjs/core";
import React from "react";

/* styles */
import styles from "./datasetSelector.css";

const buildDatasetMenuItems = (datasets) => {
  /*
  map dataset to menu item
   */
  return datasets.map((dataset) => {
    return (
      <MenuItem
        key={dataset.id}
        onClick={dataset.onClick}
        text={dataset.name}
      />
    );
  });
};
/*
 dataset menu, toggled from dataset name in app-level breadcrumbs
 */
const DatasetMenu = React.memo(({ children, datasets }) => {
  return (
    <Popover
      boundary="viewport"
      content={
        <Menu
          style={{
            maxHeight:
              "290px" /* show 9.5 datasets at 30px height each, plus top padding of 5px */,
            maxWidth: "680px" /* TODO(cc) revisit max-width versus width */,
            overflow: "auto",
          }}
        >
          {buildDatasetMenuItems(datasets)}
        </Menu>
      }
      hasBackdrop
      minimal
      modifiers={{ offset: { offset: "0, 10" } }}
      popoverClassName={styles.datasetPopover}
      position={Position.BOTTOM_LEFT}
      targetClassName={styles.datasetPopoverTarget}
    >
      {children}
    </Popover>
  );
});
export default DatasetMenu;
