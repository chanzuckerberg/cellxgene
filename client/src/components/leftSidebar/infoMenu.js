import React from "react";
import {
  Button,
  Classes,
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

const InformationMenu = React.memo((props) => {
  const {
    dataSourceLink,
    libraryVersions,
    tosURL,
    privacyURL,
    skeleton,
  } = props;
  return (
    <Popover
      content={
        <Menu>
          <MenuItem
            href="https://chanzuckerberg.github.io/cellxgene/"
            icon={<IconDocument />}
            rel="noopener"
            target="_blank"
            text="Documentation"
          />
          <MenuItem
            href="https://join-cellxgene-users.herokuapp.com/"
            icon={<IconSlack />}
            rel="noopener"
            target="_blank"
            text="Slack"
          />
          <MenuItem
            href="https://github.com/chanzuckerberg/cellxgene"
            icon={<IconGitHub />}
            rel="noopener"
            target="_blank"
            text="GitHub"
          />
          {dataSourceLink?.url && (
            <MenuItem
              href={dataSourceLink?.url}
              icon={<IconDocument />}
              rel="noopener"
              target="_blank"
              text={dataSourceLink?.name}
            />
          )}
          <MenuItem
            icon={<IconAbout />}
            popoverProps={{ openOnTargetFocus: false }}
            text="About"
          >
            <MenuItem text={libraryVersions?.cellxgene || null} />
            <MenuItem text="MIT License" />
            {tosURL && (
              <MenuItem
                href={tosURL}
                rel="noopener"
                target="_blank"
                text="Terms of Service"
              />
            )}
            {privacyURL && (
              <MenuItem
                href={privacyURL}
                rel="noopener"
                target="_blank"
                text="Privacy Policy"
              />
            )}
          </MenuItem>
        </Menu>
      }
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
        type="button"
        className={skeleton ? Classes.SKELETON : null}
      />
    </Popover>
  );
});

export default InformationMenu;
