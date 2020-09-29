import React from "react";
import { AnchorButton, ButtonGroup, Tooltip } from "@blueprintjs/core";
import * as globals from "../../globals";
import styles from "./menubar.css";

const Auth = React.memo((props) => {
  const { auth, userinfo } = props;

  if (!auth || (auth && !auth.requires_client_login)) return null;

  return (
    <ButtonGroup className={styles.menubarButton}>
      <Tooltip
        content="Log in or log out of cellxgene"
        position="bottom"
        hoverOpenDelay={globals.tooltipHoverOpenDelay}
      >
        <AnchorButton
          type="button"
          data-testid="auth-button"
          disabled={false}
          icon={!userinfo.is_authenticated ? "log-in" : "log-out"}
          href={!userinfo.is_authenticated ? auth.login : auth.logout}
        >
          {!userinfo.is_authenticated ? "Log In" : "Log Out"}
        </AnchorButton>
      </Tooltip>
    </ButtonGroup>
  );
});

export default Auth;
