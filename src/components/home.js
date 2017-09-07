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
      showClustersPlaceholderImage: false,
      heatmap: null,
      graph: null,
    };
  }
  componentWillMount() {
    // d3.csv("http://ec2-34-234-75-156.compute-1.amazonaws.com/api/v0.1/metadata", (d) => {
    //   console.log("from ec2",d)
    // })

    const prefix = "http://ec2-34-234-75-156.compute-1.amazonaws.com/api/";
    const version = "v0.1/";

    const expressions = fetch(`${prefix}${version}expression`, {
      method: "post",
      headers: new Headers({
        'Content-Type': 'application/json'
      }),
      body: JSON.stringify({
        "celllist": ["1001000173.G8", "1001000173.D4"],
        "genelist": ["1/2-SBSRNA4", "A1BG", "A1BG-AS1", "A1CF", "A2LD1", "A2M", "A2ML1", "A2MP1", "A4GALT"]
      })
    })
    // const expressions = fetch(`${prefix}${version}${expression}`)
      .then((res) => res.json())
      .then((data) => { this.setState({heatmap: data}) })

    const graph = fetch(`${prefix}${version}graph`)
      .then((res) => res.json())
      .then((data) => { this.setState({graph: data}) })

  }
  render() {
    return (
      <Container>
        <Helmet title="cellxgene" />
        <h1><Link to="/page2">cellxgene</Link></h1>
        <h3 style={{marginTop: 50}}> Gene selection criteria </h3>
        <Joy data={this.state.heatmap && this.state.heatmap.data}/>

        <Categorical/>
        <Continuous/>
        <h3 style={{marginTop: 50}}> T-SNE plot </h3>
        <button onClick={() => { this.setState({showClustersPlaceholderImage: true}) }} style={{marginBottom: 20}} className={buttonStyles.primaryButton}>Compute clustering using [n] cells in current metadata selection</button>
        <Graph data={this.state.graph}/>
        {this.state.showClustersPlaceholderImage ? <img width={750} src="https://s10.postimg.org/s6z3zbcvt/tsne.png"></img> : null}
      </Container>
    )
  }
};

export default Home;

// <Categorical title={"Sample type"} category={types}/>
// <Categorical title={"Selection"} category={selection}/>
// <Categorical title={"Location"} category={location}/>
// <Categorical title={"Sample name"} category={names}/>
