import React from "react";
import { Link } from "react-router-dom";
import _ from "lodash";
import Helmet from "react-helmet";
import Container from "./container";
import buttonStyles from "./buttons.css";

import Categorical from "./categorical/categorical";
import Continuous from "./continuous/continuous";
import Joy from "./joy/joy";
import Graph from "./graph/graph";

class Home extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      expressions: null,
      vertices: null,
      coloring: null,
      metadata: null,
    };
  }
  componentDidMount() {

    const prefix = "http://ec2-34-234-75-156.compute-1.amazonaws.com/api/";
    const version = "v0.1/";

    const expressions = fetch(`${prefix}${version}expression`, {
      method: "post",
      headers: new Headers({
        'Content-Type': 'application/json'
      }),
      body: JSON.stringify({
        // "celllist": ["1001000173.G8", "1001000173.D4"],
        "genelist": ["1/2-SBSRNA4", "A1BG", "A1BG-AS1", "A1CF", "A2LD1", "A2M", "A2ML1", "A2MP1", "A4GALT"]
      })
    })
    // const expressions = fetch(`${prefix}${version}${expression}`)
      .then((res) => res.json())
      .then((data) => { this.setState({expressions: data}) })

    const graph = fetch(`${prefix}${version}graph`)
      .then((res) => res.json())
      .then((data) => { this.setState({vertices: data}) })

    const metadata = d3.csv(`${prefix}${version}metadata`, (data) => {
      this.setState({metadata: data})
    })
  }

  createExpressionsCountsMap () {

    const CHANGE_ME_MAGIC_GENE_INDEX = 3;

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

    return (
      <Container>
        <Helmet title="cellxgene" />
        <h1><Link to="/page2">cellxgene</Link></h1>
        Showing a cluster informed downsampling of 10,000 cells.
        <button
          className={buttonStyles.primaryButton}>
           Refresh.
        </button>
        <h3 style={{marginTop: 50}}> Gene selection criteria </h3>
        {false ? <Joy data={this.state.expressions && this.state.expressions.data}/> : ""}

        <Categorical/>
        <Continuous/>
        <button
          style={{marginBottom: 20}}
          className={buttonStyles.primaryButton}>
          Compute clustering using [n] cells in current metadata selection
        </button>
        <Graph
          vertices={this.state.vertices}
          expressions={this.state.expressions}
          expressionsCountsMap={this.state.expressions && this.state.expressions ? this.createExpressionsCountsMap() : null}
          />
      </Container>
    )
  }
};

export default Home;

// <Categorical title={"Sample type"} category={types}/>
// <Categorical title={"Selection"} category={selection}/>
// <Categorical title={"Location"} category={location}/>
// <Categorical title={"Sample name"} category={names}/>
