// jshint esversion: 6
import React from "react";
import { Button, Popover, Menu, MenuItem, Position } from "@blueprintjs/core";

function InformationMenu(props) {
  const { libraryVersions, aboutLink } = props;
  return (
    <div style={{ marginLeft: 10 }} className="bp3-button-group">
      <Popover
        content={
          <Menu>
            {aboutLink ? (
              <MenuItem
                href={aboutLink}
                target="_blank"
                icon="document-open"
                text="About this dataset"
              />
            ) : (
              ""
            )}

            <MenuItem
              href="https://chanzuckerberg.github.io/cellxgene/faq.html"
              target="_blank"
              icon="help"
              text="FAQ"
            />
            <MenuItem
              href="https://join-cellxgene-users.herokuapp.com/"
              target="_blank"
              icon="chat"
              text="Chat"
            />
            <MenuItem
              href="https://chanzuckerberg.github.io/cellxgene/"
              target="_blank"
              icon="book"
              text="Docs"
            />
            <MenuItem
              href="https://github.com/chanzuckerberg/cellxgene"
              target="_blank"
              icon="git-branch"
              text="Github"
            />
            <MenuItem
              target="_blank"
              text={`cellxgene v${
                libraryVersions && libraryVersions.cellxgene
                  ? libraryVersions.cellxgene
                  : null
              }`}
            />
            <MenuItem text="MIT License" />
          </Menu>
        }
        position={Position.BOTTOM_RIGHT}
      >
        <Button
          type="button"
          className="bp3-button bp3-icon-info-sign"
          style={{
            cursor: "pointer"
          }}
        />
      </Popover>
    </div>
  );
}

export default InformationMenu;
