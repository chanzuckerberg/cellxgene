import React from "react";
import {
  Button,
  Colors,
  Menu,
  MenuItem,
  Popover,
  Position,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

/* app dependencies */
import IconAbout from "./iconAbout";
import IconDocument from "./iconDocument";
import IconGitHub from "./iconGitHub";
import IconSlack from "./iconSlack";

/* styles */
import styles from "./infoMenu.css";

const InformationMenu = React.memo((props) => {
  const { dataSourceLink, libraryVersions, tosURL, privacyURL } = props;
  return (
    <Popover
      content={
        <Menu style={{ color: Colors.BLACK }}>
          <MenuItem
            href="https://chanzuckerberg.github.io/cellxgene/"
            icon={<IconDocument />}
            rel="noopener"
            style={{ borderRadius: 0 }}
            target="_blank"
            text="Documentation"
          />
          <MenuItem
            icon={<IconSlack />}
            href="https://join-cellxgene-users.herokuapp.com/"
            rel="noopener"
            style={{ borderRadius: 0 }}
            target="_blank"
            text="Slack"
          />
          <MenuItem
            href="https://github.com/chanzuckerberg/cellxgene"
            icon={<IconGitHub />}
            rel="noopener"
            style={{ borderRadius: 0 }}
            target="_blank"
            text="GitHub"
          />
          {dataSourceLink?.url && (
            <MenuItem
              href={dataSourceLink?.url}
              icon={<IconDocument />}
              rel="noopener"
              style={{ borderRadius: 0 }}
              target="_blank"
              text={dataSourceLink?.name}
            />
          )}
          <MenuItem
            icon={<IconAbout />}
            popoverProps={{ openOnTargetFocus: false }}
            style={{
              borderRadius: 0,
              padding:
                "5px 4px 5px 7px" /* positions sub menu icon closer to menu item bounds */,
            }}
            text="About"
          >
            <MenuItem
              style={{
                borderRadius: 0,
                pointerEvents:
                  "none" /* remove hover from menu item without a link */,
              }}
              text={libraryVersions?.cellxgene || null}
            />
            <MenuItem
              style={{
                borderRadius: 0,
                pointerEvents:
                  "none" /* remove hover from menu item without a link */,
              }}
              text="MIT License"
            />
            {tosURL && (
              <MenuItem
                href={tosURL}
                rel="noopener"
                style={{ borderRadius: 0 }}
                target="_blank"
                text="Terms of Service"
              />
            )}
            {privacyURL && (
              <MenuItem
                href={privacyURL}
                rel="noopener"
                style={{ borderRadius: 0 }}
                target="_blank"
                text="Privacy Policy"
              />
            )}
          </MenuItem>
        </Menu>
      }
      popoverClassName={styles.infoMenuPopover}
      position={Position.BOTTOM_LEFT}
      modifiers={{
        hide: { enabled: false },
        preventOverflow: { enabled: false },
      }}
    >
      <Button
        data-testid="menu"
        icon={IconNames.MENU}
        minimal
        style={{
          margin: "0 2px",
          /* right padding centre aligns menu with category "color by" buttons */
          /* left padding centre aligns popover arrow with target */
        }}
        type="button"
      />
    </Popover>
  );
});

export default InformationMenu;
