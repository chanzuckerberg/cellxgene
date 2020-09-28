// jshint esversion: 6
import React from "react";
import { Button, Popover, Menu, MenuItem, Position } from "@blueprintjs/core";
import styles from "./menubar.css";

const handleClick = (dispatch) => {
  dispatch({ type: "toggle dataset drawer" });
};

const InformationMenu = React.memo((props) => {
  const {
    libraryVersions,
    tosURL,
    privacyURL,
    auth,
    userinfo,
    dispatch,
  } = props;
  return (
    <div className={`bp3-button-group ${styles.menubarButton}`}>
      <Popover
        content={
          <Menu>
            <MenuItem
              onClick={() => handleClick(dispatch)}
              icon="info-sign"
              text="Dataset Overview"
            />

            <MenuItem
              href="https://chanzuckerberg.github.io/cellxgene/"
              target="_blank"
              icon="book"
              text="Documentation"
            />
            <MenuItem
              href="https://join-cellxgene-users.herokuapp.com/"
              target="_blank"
              icon="chat"
              text="Chat"
            />
            <MenuItem
              href="https://github.com/chanzuckerberg/cellxgene"
              target="_blank"
              icon="git-branch"
              text="Github"
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
              />
            ) : null}

            {auth.requires_client_login && userinfo.is_authenticated ? (
              <MenuItem text={`Logged in as: ${userinfo.email}`} />
            ) : null}
            {auth.requires_client_login && userinfo.is_authenticated ? (
              <MenuItem text="Log Out" href={auth.logout} />
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
          type="button"
          className="bp3-button bp3-icon-info-sign"
          style={{
            cursor: "pointer",
          }}
        />
      </Popover>
    </div>
  );
});

export default InformationMenu;
