import React from "react";
import {
  AnchorButton,
  Button,
  MenuItem,
  Tooltip,
  Popover,
  Menu,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

import * as globals from "../../globals";
import styles from "./menubar.css";

const BASE_EMOJI = [0x1f9d1, 0x1f468, 0x1f469];
const SKIN_TONES = [0x1f3fb, 0x1f3fc, 0x1f3fd, 0x1f3fe, 0x1f3ff];
const MICROSCOPE = 0x1f52c;
const ZERO_WIDTH_JOINER = 0x0200d;

const Auth = React.memo((props) => {
  const { auth, userinfo } = props;

  const randomInt = Math.random() * 15;
  const sexIndex = Math.floor(randomInt / 5);
  const skinToneIndex = Math.floor(randomInt % 5);

  const scientist = String.fromCodePoint(
    BASE_EMOJI[sexIndex],
    SKIN_TONES[skinToneIndex],
    ZERO_WIDTH_JOINER,
    MICROSCOPE
  );

  if (!auth?.["requires_client_login"]) return null;

  if (userinfo?.["is_authenticated"]) {
    const PopoverContent = (
      <Menu>
        <MenuItem
          data-testid="user-email"
          text={`Logged in as: ${userinfo.email}`}
        />
        <MenuItem
          data-testid="log-out"
          text="Log Out"
          href={auth.logout}
          icon={IconNames.LOG_OUT}
        />
      </Menu>
    );

    return (
      <Popover content={PopoverContent}>
        <Button
          data-testid="user-info"
          className={styles.menubarButton}
          style={{ padding: 0 }}
        >
          {/*  eslint-disable-next-line no-constant-condition -- disable profile picture until CSP is tweaked */}
          {userinfo?.picture && false ? (
            <img alt="profile" size="21px" src={userinfo?.picture} />
          ) : (
            <span style={{ fontSize: "18px" }}>{scientist}</span>
          )}
        </Button>
      </Popover>
    );
  }

  return (
    <Tooltip
      content="Log in to cellxgene"
      position="bottom"
      hoverOpenDelay={globals.tooltipHoverOpenDelay}
    >
      <AnchorButton
        type="button"
        data-testid="log-in"
        disabled={false}
        href={auth.login}
        className={styles.menubarButton}
      >
        Log In
      </AnchorButton>
    </Tooltip>
  );
});

export default Auth;
