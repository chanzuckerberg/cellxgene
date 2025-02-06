
import React from "react";
import { connect } from "react-redux";

import { FilterOutlined, HistoryOutlined } from "@ant-design/icons";
import { Button } from "antd";
import * as globals from "../../globals";
import cls from "./layout.css";

@connect((state) => ({
  pin: state.pin,
}))
class Layout extends React.Component {
  /*
    Layout - this react component contains all the layout style and logic for the application once it has loaded.

    The layout is based on CSS grid: the left and right sidebars have fixed widths, the graph in the middle takes the
    remaining space.

    Note, the renderGraph child is a function rather than a fully-instantiated element because the middle pane of the
    app is dynamically-sized. It must have access to the containing viewport in order to know how large the graph
    should be.
  */
  componentDidMount() {
    /*
      This is a bit of a hack. In order for the graph to size correctly, it needs to know the size of the parent
      viewport. Unfortunately, it can only do this once the parent div has been rendered, so we need to render twice.
    */
    this.forceUpdate();
  }

  render() {
    const { children, pin, dispatch } = this.props;
    const [leftSidebar, renderGraph, rightSidebar] = children;

    return (
      <div className={cls.root}>
        <div className={`${cls.left} ${pin.left ? cls.pinned : ""}`}>
          {leftSidebar}
          {!pin.left && (
            <Button
              className={cls.trigger}
              type="dashed"
              icon={<FilterOutlined />}
              onClick={() => {
                dispatch({
                  type: "pin: update",
                  loc: "left",
                  pinned: true,
                });
              }}
            />
          )}
        </div>
        <div
          className={cls.graph}
          ref={(ref) => {
            this.viewportRef = ref;
          }}
        >
          {this.viewportRef ? renderGraph(this.viewportRef) : null}
        </div>
        <div
          className={`${cls.right} ${pin.right ? cls.pinned : ""}`}
          style={{ borderLeft: `1px solid ${globals.lightGrey}` }}
        >
          {rightSidebar}
          {!pin.right && (
            <Button
              className={cls.trigger}
              type="dashed"
              icon={<HistoryOutlined />}
              onClick={() => {
                dispatch({
                  type: "pin: update",
                  loc: "right",
                  pinned: true,
                });
              }}
            />
          )}
        </div>
      </div>
    );
  }
}

export default Layout;
