// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import HistogramBrush from "../brushableHistogram";
import { Button, ButtonGroup, Tooltip } from "@blueprintjs/core";
import actions from "../../actions";

@connect((state) => {
  return {};
})
class Gene extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const { gene } = this.props;
    return (
      <div>
        <p> {gene} </p>
      </div>
    );
  }
}

export default Gene;
