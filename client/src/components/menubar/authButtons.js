import React from "react";
import { AnchorButton, Tooltip } from "@blueprintjs/core";
import * as globals from "../../globals";
import styles from "./menubar.css";

const Auth = React.memo((props) => {
  const { auth } = props;

  if (!auth || (auth && !auth.requires_client_login)) return null;

  return (
    <div className={`bp3-button-group ${styles.menubarButton}`}>
      <Tooltip
        content="Log in or log out of cellxgene"
        position="bottom"
        hoverOpenDelay={globals.tooltipHoverOpenDelay}
      >
        <AnchorButton
          type="button"
          data-testid="auth-button"
          disabled={false}
          icon={!auth.is_authenticated ? "log-in" : "log-out"}
          href={!auth.is_authenticated ? auth.login : auth.logout}
        >
          {!auth.is_authenticated ? "Log In" : "Log Out"}
        </AnchorButton>
      </Tooltip>
    </div>
  );
});

export default Auth;
