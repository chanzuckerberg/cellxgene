import React from "react";
import { Button, Menu, MenuItem, Popover, Position } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

/* app dependencies */
import IconAbout from "./iconAbout";
import IconDocument from "./iconDocument";
import IconGitHub from "./iconGitHub";
import IconSlack from "./iconSlack";

const InformationMenu = React.memo((props) => {
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'libraryVersions' does not exist on type ... Remove this comment to see the full error message
  const { libraryVersions, tosURL, privacyURL } = props;
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
            target="_blank"
            text="Chat"
            rel="noopener"
          />
          <MenuItem
            href="https://github.com/chanzuckerberg/cellxgene"
            icon={<IconGitHub />}
            target="_blank"
            text="Github"
            rel="noopener"
          />
          <MenuItem
            icon={<IconAbout />}
            popoverProps={{ openOnTargetFocus: false }}
            text="About cellxgene"
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
        preventOverflow: { enabled: false },
        hide: { enabled: false },
      }}
    >
      <Button data-testid="menu" icon={IconNames.MENU} minimal type="button" />
    </Popover>
  );
});

export default InformationMenu;
