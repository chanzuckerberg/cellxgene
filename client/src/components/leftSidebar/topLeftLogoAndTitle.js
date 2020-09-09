// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

import * as globals from "../../globals";
import Logo from "../framework/logo";
import Truncate from "../util/truncate";
import InfoDrawer from "../infoDrawer/infoDrawer";

const DATASET_TITLE_FONT_SIZE = 14;

@connect((state) => ({
  datasetTitle: state.config?.displayNames?.dataset ?? "",
  hoverState: state.controls.singletonHover,
}))
class LeftSideBar extends React.Component {
  handleClick = () => {
    const { dispatch } = this.props;
    dispatch({ type: "toggle dataset drawer" });
  };

  render() {
    const { datasetTitle, hoverState } = this.props;

    return (
      <div
        style={{
          paddingLeft: 8,
          paddingTop: 8,
          width: globals.leftSidebarWidth,
          zIndex: 1,
          borderBottom: `1px solid ${globals.lighterGrey}`,
        }}
      >
        <Logo size={30} />
        <span
          style={{
            fontSize: 28,
            position: "relative",
            top: -6,
            fontWeight: "bold",
            marginLeft: 5,
            color: globals.logoColor,
            userSelect: "none",
          }}
        >
          cell
          <span
            style={{
              position: "relative",
              top: 1,
              fontWeight: 300,
              fontSize: 24,
            }}
          >
            ×
          </span>
          gene
        </span>
        <Button
          outlined
          icon={IconNames.BOOK}
          style={{
            fontSize: DATASET_TITLE_FONT_SIZE,
            marginLeft: "7px",
            position: "relative",
            top: -8,
          }}
          active={hoverState}
          onClick={this.handleClick}
        >
          <Truncate>
            <span style={{ maxWidth: 155 }} data-testid="header">
              {datasetTitle}
            </span>
          </Truncate>
        </Button>
        <InfoDrawer />
      </div>
    );
  }
}

export default LeftSideBar;
