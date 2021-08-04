import React from "react";
import Helmet from "react-helmet";
import { connect } from "react-redux";

import Container from "./framework/container";
import Layout from "./framework/layout";
import LeftSideBar from "./leftSidebar";
import RightSideBar from "./rightSidebar";
import Legend from "./continuousLegend";
import Graph from "./graph/graph";
import Dotplot from "./dotplot";
import MenuBar from "./menubar";
import Autosave from "./autosave";
import Embedding from "./embedding";
import TermsOfServicePrompt from "./termsPrompt";

import actions from "../actions";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  loading: (state as any).controls.loading,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  error: (state as any).controls.error,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  graphRenderCounter: (state as any).controls.graphRenderCounter,
  layoutChoice: (state as any).layoutChoice,
}))
class App extends React.Component {
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  componentDidMount() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    /* listen for url changes, fire one when we start the app up */
    window.addEventListener("popstate", this._onURLChanged);
    this._onURLChanged();
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 0 arguments, but got 1.
    dispatch(actions.doInitialDataLoad(window.location.search));
    this.forceUpdate();
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  _onURLChanged() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch({ type: "url changed", url: document.location.href });
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'loading' does not exist on type 'Readonl... Remove this comment to see the full error message
    const { loading, error, graphRenderCounter, layoutChoice } = this.props;
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
            error loading cellxgene
          </div>
        ) : null}
        {loading || error ? null : (
          <Layout>
            <LeftSideBar />
            {/* eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS. */}
            {(viewportRef: any) => (
              <>
                <MenuBar />
                <Embedding />
                <Autosave />
                <TermsOfServicePrompt />
                {/* @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call. */}
                <Legend viewportRef={viewportRef} />
                {layoutChoice.dotplot && <Dotplot viewportRef={viewportRef} />}

                {/* @ts-expect-error ts-migrate(2322) FIXME: Type '{ key: any; viewportRef: any; }' is not assi... Remove this comment to see the full error message */}
                <Graph
                  key={graphRenderCounter}
                  dotplotMode={layoutChoice.dotplot}
                  viewportRef={viewportRef}
                />
              </>
            )}
            <RightSideBar />
          </Layout>
        )}
      </Container>
    );
  }
}

export default App;
