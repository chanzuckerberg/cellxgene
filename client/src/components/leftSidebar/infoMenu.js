// jshint esversion: 6
import React from "react";
import { Button, Menu, MenuItem, Popover, Position } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

const InformationMenu = React.memo((props) => {
  const { libraryVersions, tosURL, privacyURL } = props;
  return (
    <Popover
      content={
        <Menu>
          <MenuItem
            href="https://chanzuckerberg.github.io/cellxgene/"
            target="_blank"
            icon="book"
            text="Documentation"
            rel="noopener"
          />
          <MenuItem
            href="https://join-cellxgene-users.herokuapp.com/"
            target="_blank"
            icon="chat"
            text="Chat"
            rel="noopener"
          />
          <MenuItem
            href="https://github.com/chanzuckerberg/cellxgene"
            target="_blank"
            icon="git-branch"
            text="Github"
            rel="noopener"
          />
          <MenuItem
            target="_blank"
            text={
              libraryVersions && libraryVersions.cellxgene
                ? libraryVersions.cellxgene
                : null
            }
          />
          <MenuItem text="MIT License" />
          {tosURL ? (
            <MenuItem href={tosURL} target="_blank" text="Terms of Service" />
          ) : null}
          {privacyURL ? (
            <MenuItem
              href={privacyURL}
              target="_blank"
              text="Privacy Policy"
              rel="noopener"
            />
          ) : null}
        </Menu>
      }
      position={Position.BOTTOM_RIGHT}
      modifiers={{
        preventOverflow: { enabled: false },
        hide: { enabled: false },
      }}
    >
      <Button
        data-testid="menu"
        type="button"
        icon={IconNames.INFO_SIGN}
        style={{
          cursor: "pointer",
          verticalAlign: "middle",
        }}
      />
    </Popover>
  );
});

export default InformationMenu;
