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

    /* listen for resize events */
    window.addEventListener("resize", () => {
      dispatch({
        type: "window resize",
        data: {
          height: window.innerHeight,
          width: window.innerWidth,
        },
      });
    });
    dispatch({
      type: "window resize",
      data: {
        height: window.innerHeight,
        width: window.innerWidth,
      },
    });
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
        <div>
          {loading ? null : <LeftSideBar />}
          {loading ? null : <RightSideBar />}
          {loading ? null : <MenuBar />}
          {loading ? null : <Graph key={graphRenderCounter} />}
          {loading ? null : <Autosave />}
          {loading ? null : <TermsOfServicePrompt />}
          <Legend />
        </div>
      </Container>
    );
  }
}

export default App;
