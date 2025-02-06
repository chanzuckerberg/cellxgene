import React from "react";
import Helmet from "react-helmet";
import { connect } from "react-redux";

import Container from "./framework/container";
import Layout from "./framework/layout";
import LeftSideBar from "./leftSidebar";
import Legend from "./continuousLegend";
import Graph from "./graph/graph";
import MenuBar from "./menubar";
// import Autosave from "./autosave";
import Embedding from "./embedding";
import Conversation from "./conversation";

import actions from "../actions";

@connect((state) => ({
  loading: state.controls.loading,
  error: state.controls.error,
  graphRenderCounter: state.controls.graphRenderCounter,
}))
class App extends React.Component {
  componentDidMount() {
    const { dispatch } = this.props;

    /* listen for url changes, fire one when we start the app up */
    window.addEventListener("popstate", this._onURLChanged);
    this._onURLChanged();

    dispatch(actions.doInitialDataLoad(window.location.search));
    this.forceUpdate();
  }

  _onURLChanged() {
    const { dispatch } = this.props;

    dispatch({ type: "url changed", url: document.location.href });
  }

  render() {
    const { loading, error, graphRenderCounter } = this.props;
    return (
      <Container>
        <Helmet title="Dandelion" />
        {loading ? (
          <div
            style={{
              position: "fixed",
              fontWeight: 500,
              top: window.innerHeight / 2,
              left: window.innerWidth / 2 - 50,
            }}
          >
            loading
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
        {loading || error ? null : (
          <Layout>
            <LeftSideBar />
            {(viewportRef) => (
              <>
                <MenuBar />
                <Embedding />
                {/* <Autosave /> */}
                <Legend viewportRef={viewportRef} />
                <Graph key={graphRenderCounter} viewportRef={viewportRef} />
              </>
            )}
            <Conversation />
          </Layout>
        )}
      </Container>
    );
  }
}

export default App;
