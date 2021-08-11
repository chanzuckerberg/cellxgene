/* core dependencies */
import { Menu, MenuItem, Popover, Position } from "@blueprintjs/core";
import React from "react";

/* styles */
// @ts-expect-error --- TODO fix import
import styles from "./datasetSelector.css";

// @ts-expect-error --- TODO add typing for datasets
const buildDatasetMenuItems = (datasets) =>
  /*
    map dataset to menu item
     */
  // @ts-expect-error --- TODO add typing for dataset
  datasets.map((dataset) => (
    <MenuItem key={dataset.id} onClick={dataset.onClick} text={dataset.name} />
  ));

/*
 dataset menu, toggled from dataset name in app-level breadcrumbs
 */
// @ts-expect-error --- TODO add typing for props
const DatasetMenu = React.memo(({ children, datasets }) => (
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
));
export default DatasetMenu;
