// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import { useAuth0 } from "../../react-auth0-spa";
import Categorical from "../categorical";
import * as globals from "../../globals";
import DynamicScatterplot from "../scatterplot/scatterplot";
import TopLeftLogoAndTitle from "./topLeftLogoAndTitle";
import AuthN from "./AuthN";

const LeftSideBar = (props) => {
  const { scatterplotXXaccessor, scatterplotYYaccessor } = props;
  const { loading } = useAuth0();

  if (loading) {
    return <div>Loading...</div>;
  }

  return (
    <div
      style={{
        /* x y blur spread color */
        borderRight: `1px solid ${globals.lightGrey}`,
        display: "flex",
        flexDirection: "column",
        height: "100%",
      }}
    >
      <TopLeftLogoAndTitle />
      <div
        style={{
          height: "100%",
          width: globals.leftSidebarWidth,
          overflowY: "auto",
        }}
      >
        <AuthN />
        <Categorical />
      </div>
      {scatterplotXXaccessor && scatterplotYYaccessor ? (
        <DynamicScatterplot />
      ) : null}
    </div>
  );
};

const mapStateToProps = (state) => ({
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
});

export default connect(mapStateToProps)(LeftSideBar);
