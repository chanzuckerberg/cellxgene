import React from "react";
import _ from "lodash";
import Categorical from "./categorical/categorical";
import Expression from "./expression/expression";
import { connect } from "react-redux";
import Heatmap from "./expression/diffExpHeatmap";
import * as globals from "../globals";

@connect((state) => {
  return {

  }
})
class LeftSideBar extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      currentTab: "metadata"
    };
  }
  render() {
    return (
      <div style={{position: "fixed"}}>
        <div style={{padding: 10}}>
          <button
            style={{
              padding: "10px 30px",
              outline: 0,
              fontSize: 18,
              fontStyle: this.state.currentTab === "metadata" ? "inherit" : "italic",
              cursor: "pointer",
              border: "none",
              backgroundColor: "#FFF",
              borderTop: this.state.currentTab === "metadata" ? "4px solid " + globals.brightBlue : "none",
              borderBottom: this.state.currentTab === "metadata" ? "none" : "1px solid " + globals.lightGrey,
              borderRight: this.state.currentTab === "metadata" ? "1px solid " + globals.lightGrey : "none",
              borderLeft: this.state.currentTab === "metadata" ? "1px solid " + globals.lightGrey : "none",
            }}
            onClick={() => {this.setState({currentTab: "metadata"})}}>
            Metadata
          </button>
          <button
            style={{
              padding: "10px 30px",
              outline: 0,
              fontSize: 18,
              fontStyle: this.state.currentTab === "expression" ? "inherit" : "italic",
              cursor: "pointer",
              border: "none",
              backgroundColor: "#FFF",
              borderTop: this.state.currentTab === "expression" ? "4px solid " + globals.brightBlue : "none",
              borderBottom: this.state.currentTab === "expression" ? "none" : "1px solid " + globals.lightGrey,
              borderRight: this.state.currentTab === "expression" ? "1px solid " + globals.lightGrey : "none",
              borderLeft: this.state.currentTab === "expression" ? "1px solid " + globals.lightGrey : "none",
            }}
            onClick={() => {this.setState({currentTab: "expression"})}}>
            Expression
          </button>
        </div>
        <div style={{
          height: 500,
          width: 350,
          padding: 10,
          overflowY: "scroll",
          overflowX: "hidden",
        }}>
          {this.state.currentTab === "metadata" ? <Categorical/> : null}
          {this.state.currentTab === "expression" ? <Heatmap/> : null}
        </div>
        <div style={{
          boxShadow: "-3px -4px 13px 0px rgba(201,201,201,1)",
          paddingTop: 10
        }}>
          <Expression/>
        </div>
      </div>
    )
  }
};

export default LeftSideBar;
