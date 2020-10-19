import React, { useState } from "react";

import {
  AnchorButton,
  Button,
  MenuItem,
  Tooltip,
  Popover,
  Menu,
  Elevation,
  PopoverPosition,
  Checkbox,
  Card,
} from "@blueprintjs/core";

import { IconNames } from "@blueprintjs/icons";

import * as globals from "../../globals";

import styles from "./menubar.css";

import { storageGet, storageSet, KEYS } from "../util/localStorage";

const BASE_EMOJI = [0x1f9d1, 0x1f468, 0x1f469];
const SKIN_TONES = [0x1f3fb, 0x1f3fc, 0x1f3fd, 0x1f3fe, 0x1f3ff];
const MICROSCOPE = 0x1f52c;
const ZERO_WIDTH_JOINER = 0x0200d;

const LOGIN_PROMPT_OFF = "off";

const Auth = React.memo((props) => {
  const [isPromptOpen, setIsPromptOpen] = useState(shouldShowPrompt());

  const { auth, userInfo } = props;

  const isAuthenticated = userInfo && userInfo.is_authenticated;

  window.userInfo = userInfo;

  const randomInt = Math.random() * 15;
  const sexIndex = Math.floor(randomInt / 5);
  const skinToneIndex = Math.floor(randomInt % 5);

  const scientist = String.fromCodePoint(
    BASE_EMOJI[sexIndex],
    SKIN_TONES[skinToneIndex],
    ZERO_WIDTH_JOINER,
    MICROSCOPE
  );

  if (!shouldShowAuth()) return null;

  if (isAuthenticated) {
    const PopoverContent = (
      <Menu>
        <MenuItem
          data-testid="user-email"
          text={`Logged in as: ${userInfo.email}`}
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
          {userInfo?.picture && false ? (
            <img alt="profile" size="21px" src={userInfo?.picture} />
          ) : (
            <span style={{ fontSize: "18px" }}>{scientist}</span>
          )}
        </Button>
      </Popover>
    );
  }

  const LoginButton = (
    <Tooltip
      content="Log in to cellxgene"
      position="bottom"
      hoverOpenDelay={globals.tooltipHoverOpenDelay}
    >
      <AnchorButton
        type="button"
        data-testid="log-in"
        href={auth.login}
        className={styles.menubarButton}
      >
        Log In
      </AnchorButton>
    </Tooltip>
  );

  if (isPromptOpen) {
    return (
      <Popover
        position={PopoverPosition.AUTO_END}
        isOpen
        content={<PromptContent setIsPromptOpen={setIsPromptOpen} />}
        onInteraction={setIsPromptOpen}
      >
        {LoginButton}
      </Popover>
    );
  }

  return LoginButton;

  function shouldShowAuth() {
    return auth && auth.requires_client_login;
  }

  function shouldShowPrompt() {
    if (storageGet(KEYS.LOGIN_PROMPT) === LOGIN_PROMPT_OFF) return false;

    return shouldShowAuth && !isAuthenticated;
  }
});

function PromptContent({ setIsPromptOpen }) {
  const [isChecked, setIsChecked] = useState(false);

  function handleOKClick() {
    if (isChecked) {
      storageSet(KEYS.LOGIN_PROMPT, LOGIN_PROMPT_OFF);
    }

    setIsPromptOpen(false);
  }

  function handleCheckboxChange() {
    setIsChecked(!isChecked);
  }

  return (
    <Card style={{ width: "500px" }} elevation={Elevation.TWO}>
      <p>
        Logging in will enable you to create your own categories and labels.
        Logging in later will reset cellxgene to the default view and cause you
        to lose progress.
      </p>
      <Checkbox
        style={{ width: "230px" }}
        checked={isChecked}
        onChange={handleCheckboxChange}
        data-testid="login-hint-do-not-show-again"
      >
        Do not show me this message again
      </Checkbox>
      <div
        style={{ display: "flex", justifyContent: "flex-end", marginTop: 15 }}
      >
        <Button
          onClick={handleOKClick}
          intent="primary"
          data-testid="login-hint-yes"
        >
          Acknowledge
        </Button>
      </div>
    </Card>
  );
}

export default Auth;
