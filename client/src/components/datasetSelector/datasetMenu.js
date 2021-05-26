/* core dependencies */
import { Colors, Menu, MenuItem, Popover, Position } from "@blueprintjs/core";
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
const DatasetMenu = React.memo(({ children, datasets }) => {
  const menuScrollable = datasets.length > 9; /* scrollable at 10 datasets */
  return (
    <Popover
      boundary="viewport"
      content={
        <Menu
          className={menuScrollable ? styles.datasetMenuScrollable : null}
          style={{
            color: Colors.BLACK,
            maxHeight:
              "290px" /* show 9.5 datasets at 30px height each, plus top padding of 5px */,
            maxWidth: "680px" /* TODO(cc) revisit max-width versus width */,
            overflow: "auto",
            paddingRight: menuScrollable
              ? 0 /* scrollbar adds padding to menu  */ /* TODO(cc) address ff requiring additional padding */
              : null,
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
