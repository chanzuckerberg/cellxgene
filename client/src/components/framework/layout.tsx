import React from "react";
import * as globals from "../../globals";

class Layout extends React.Component {
  viewportRef: any;
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
    const { children } = this.props;
    // @ts-expect-error ts-migrate(2488) FIXME: Type 'ReactNode' must have a '[Symbol.iterator]()'... Remove this comment to see the full error message
    const [leftSidebar, renderGraph, rightSidebar] = children;
    return (
      <div
        style={{
          display: "grid",
          gridTemplateColumns: `
          [left-sidebar-start] ${globals.leftSidebarWidth + 1}px
          [left-sidebar-end graph-start] auto
          [graph-end right-sidebar-start] ${
            globals.rightSidebarWidth + 1
          }px [right-sidebar-end]
        `,
          gridTemplateRows: "[top] auto [bottom]",
          gridTemplateAreas: "left-sidebar | graph | right-sidebar",
          columnGap: "0px",
          justifyItems: "stretch",
          alignItems: "stretch",
          height: "inherit",
          width: "inherit",
          position: "relative",
          top: 0,
          left: 0,
          minWidth: "1240px",
        }}
      >
        <div
          style={{
            gridArea: "top / left-sidebar-start / bottom / left-sidebar-end",
            position: "relative",
            height: "inherit",
            overflowY: "auto",
          }}
        >
          {leftSidebar}
        </div>
        <div
          style={{
            zIndex: 0,
            gridArea: "top / graph-start / bottom / graph-end",
            position: "relative",
            height: "inherit",
          }}
          ref={(ref) => {
            this.viewportRef = ref;
          }}
        >
          {this.viewportRef ? renderGraph(this.viewportRef) : null}
        </div>
        <div
          style={{
            gridArea: "top / right-sidebar-start / bottom / right-sidebar-end",
            position: "relative",
            height: "inherit",
            overflowY: "auto",
          }}
        >
          {rightSidebar}
        </div>
      </div>
    );
  }
}

export default Layout;
