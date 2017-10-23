import React from "react";
import _ from "lodash";
import Helmet from "react-helmet";
import Container from "./framework/container";
import buttonStyles from "./framework/buttons.css";
import { connect } from "react-redux";
import PulseLoader from "halogen/PulseLoader";

import Categorical from "./categorical/categorical";
import Continuous from "./continuous/continuous";
import Joy from "./joy/joy";
import Graph from "./graph/graph";
import * as globals from "../globals";
import * as actions from "../actions";

import SectionHeader from "./framework/sectionHeader";

@connect((state) => {
  return {
    cells: state.cells
  }
})
class App extends React.Component {
  constructor(props) {
    super(props);
    this.state = {

    };
  }
  _onURLChanged () {
    this.props.dispatch({type: 'url changed', url: document.location.href});
  };
  componentDidMount() {

    /* listen for url changes, fire one when we start the app up */
    window.addEventListener('popstate', this._onURLChanged);
    this._onURLChanged();

    this.props.dispatch(actions.requestCells())

  }

  createExpressionsCountsMap () {

    const CHANGE_ME_MAGIC_GENE_INDEX = 5;

    const expressionsCountsMap = {};

    /* currently selected gene */
    expressionsCountsMap.geneName = this.state.expressions.data.genes[3];

    let maxExpressionValue = 0;

    /* create map of expressions for every cell */
    this.state.expressions.data.cells.map((c) => {
      /* cellname = 234 */
      expressionsCountsMap[c.cellname] = c["e"][CHANGE_ME_MAGIC_GENE_INDEX];
      /* collect the maximum value as we iterate */
      if (c["e"][CHANGE_ME_MAGIC_GENE_INDEX] > maxExpressionValue) {
        maxExpressionValue = c["e"][CHANGE_ME_MAGIC_GENE_INDEX]
      }
    })

    expressionsCountsMap.maxValue = maxExpressionValue;

    return expressionsCountsMap;
  }

  render() {
    // console.log('app:', this.props, this.state)

    return (
      <Container>
        <Helmet title="cellxgene" />
        {
          this.props.cells.loading ?
          <div style={{display: "flex", justifyContent: "center", marginTop: 20}}>
            <PulseLoader color="rgb(0,0,0)" size="10px" margin="4px"/>
            <span style={{marginLeft: 10, fontFamily: globals.accentFont, fontStyle: "italic"}}>loading cells</span>
          </div> :
          null
        }
        {this.props.cells.error ? "Error loading cells" : null}
        <SectionHeader text="Gene Selection Criteria"/>
        {false ? <Joy data={this.state.expressions && this.state.expressions.data}/> : ""}
        <Categorical/>
        <Continuous/>
        <Graph
          vertices={this.state.vertices}
          expressions={this.state.expressions}
          expressionsCountsMap={this.state.expressions && this.state.expressions ? this.createExpressionsCountsMap() : null}
          />
      </Container>
    )
  }
};

export default App;
