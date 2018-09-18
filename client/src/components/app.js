// jshint esversion: 6
import React from "react";
import _ from "lodash";
import Helmet from "react-helmet";
import Container from "./framework/container";
import { connect } from "react-redux";
// import PulseLoader from "halogen/PulseLoader";

import LeftSideBar from "./leftsidebar";
import Legend from "./continuousLegend";
import Graph from "./graph/graph";
import * as globals from "../globals";
import actions from "../actions";

import SectionHeader from "./framework/sectionHeader";

@connect(state => {
  return {
    loading: state.controls.loading,
    error: state.controls.error
  };
})
class App extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  _onURLChanged() {
    this.props.dispatch({ type: "url changed", url: document.location.href });
  }

  componentDidMount() {
    /* listen for url changes, fire one when we start the app up */
    window.addEventListener("popstate", this._onURLChanged);
    this._onURLChanged();

    this.props.dispatch(actions.doInitialDataLoad(window.location.search));

    /* listen for resize events */
    window.addEventListener("resize", () => {
      this.props.dispatch({
        type: "window resize",
        data: {
          height: window.innerHeight,
          width: window.innerWidth
        }
      });
    });
    this.props.dispatch({
      type: "window resize",
      data: {
        height: window.innerHeight,
        width: window.innerWidth
      }
    });
  }

  render() {
    const { loading, error } = this.props;
    return (
      <Container>
        <Helmet title="cellxgene" />
        {loading ? (
          <div
            style={{
              position: "fixed",
              fontWeight: 500,
              top: window.innerHeight / 2,
              left: window.innerWidth / 2 - 50
            }}
          >
            loading cellxgene
          </div>
        ) : null}
        {error ? "Error loading cells" : null}
        <div>
          {loading ? null : <LeftSideBar />}
          <div
            style={{
              padding: 15,
              width: 1440 - 410 /* but responsive */,
              marginLeft: 350 /* but responsive */
            }}
          >
            {loading ? null : <Graph />}

            <Legend />
            {}
          </div>
        </div>
      </Container>
    );
  }
}

export default App;
