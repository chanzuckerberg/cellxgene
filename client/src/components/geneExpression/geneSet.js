// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import HistogramBrush from "../brushableHistogram";
import * as globals from "../../globals";
import actions from "../../actions";

@connect(state => {
  return {};
})
class GeneSet extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const { setName, setGenes } = this.props;
    return (
      <div>
        {setName}
        {_.map(setGenes, gene => {
          return <p key={gene}> {gene} </p>;
        })}
      </div>
    );
  }
}

export default GeneSet;
