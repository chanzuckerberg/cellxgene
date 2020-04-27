// jshint esversion: 6
import React from "react";
import Helmet from "react-helmet";
import { connect } from "react-redux";

import Container from "./framework/container";
import LeftSideBar from "./leftSidebar";
import RightSideBar from "./rightSidebar";
import Legend from "./continuousLegend";
import Graph from "./graph/graph";
import MenuBar from "./menubar";
import Autosave from "./autosave";
import TermsOfServicePrompt from "./termsPrompt";
import * as globals from "../globals";

import actions from "../actions";

@connect((state) => ({
  loading: state.controls.loading,
  error: state.controls.error,
  graphRenderCounter: state.controls.graphRenderCounter,
}))
class App extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  componentDidMount() {
    const { dispatch } = this.props;

    /* listen for url changes, fire one when we start the app up */
    window.addEventListener("popstate", this._onURLChanged);
    this._onURLChanged();

    dispatch(actions.doInitialDataLoad(window.location.search));
  }

  componentDidUpdate() {
    const { secondRender } = this.state;
    if (!secondRender) this.setState({ secondRender: true })
  }

  _onURLChanged() {
    const { dispatch } = this.props;

    dispatch({ type: "url changed", url: document.location.href });
  }

  render() {
    const { loading, error, graphRenderCounter } = this.props;
    return (
      <Container>
        <Helmet title="cellxgene" />
        {loading ? (
          <div
            style={{
              position: "fixed",
              fontWeight: 500,
              top: window.innerHeight / 2,
              left: window.innerWidth / 2 - 50,
            }}
          >
            loading cellxgene
          </div>
        ) : null}
        {error ? (
          <div
            style={{
              position: "fixed",
              fontWeight: 500,
              top: window.innerHeight / 2,
              left: window.innerWidth / 2 - 50,
            }}
          >
            error loading
          </div>
        ) : null}
        {loading ? null : <div
          style={{
            display: "grid",
            gridTemplateColumns: `
              [left-sidebar-start] ${globals.leftSidebarWidth+1}px
              [left-sidebar-end graph-start] auto
              [graph-end right-sidebar-start] ${globals.rightSidebarWidth+1}px [right-sidebar-end]
            `,
            gridTemplateRows: "[top] auto [bottom]",
            gridTemplateAreas: "left-sidebar | graph | right-sidebar",
            columnGap: "0px",
            justifyItems: "stretch",
            alignItems: "stretch",
            height: "100vh",
            width: "100vw",
            position: "relative",
            top: 0,
            left: 0,
            minHeight: "inherit",
            minWidth: "1240px",
          }}
        >
          <LeftSideBar
            style={{
              gridArea: "top / left-sidebar-start / bottom / left-sidebar-end",
              position: "relative"
            }}
          />
          <div
            style={{
              zIndex: 0,
              gridArea: "top / graph-start / bottom / graph-end",
              height: "100%",
              width: "100%",
              position: "relative"
            }}
            ref={ref => { this.viewportRef = ref; }}
          >
            {
              this.viewportRef
                ? <>
                  <MenuBar />
                  <Autosave />
                  <TermsOfServicePrompt />
                  <Legend viewportRef={this.viewportRef} />
                  <Graph
                    key={ graphRenderCounter }
                    viewportRef={this.viewportRef}
                  />
                </>
                : null
            }
          </div>
          <RightSideBar
            style={{
              gridArea: "top / right-sidebar-start / bottom / right-sidebar-end",
            }}
          />
        </div>}
      </Container>
    );
  }
}

export default App;
