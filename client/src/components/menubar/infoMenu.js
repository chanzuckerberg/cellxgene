// jshint esversion: 6
import React from "react";
import {
  Button,
  ButtonGroup,
  Classes,
  Menu,
  MenuItem,
  Popover,
  Position,
} from "@blueprintjs/core";
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
    <ButtonGroup className={`${styles.menubarButton}`}>
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

            {auth?.["requires_client_login"] &&
            userinfo?.["is_authenticated"] ? (
              <>
                <MenuItem text={`Logged in as: ${userinfo.email}`} />
                <MenuItem text="Log Out" href={auth.logout} />
              </>
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
          className={`${Classes.BUTTON} bp3-icon-info-sign`}
          style={{
            cursor: "pointer",
          }}
        />
      </Popover>
    </ButtonGroup>
  );
});

export default InformationMenu;
