import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import styles from "./expression.css";
import SectionHeader from "../framework/sectionHeader"

@connect((state) => {

  return {
    expression: state.expression.data
  }
})
class Expression extends React.Component {
  constructor(props) {
    super(props);
    this.state = {

    };
  }
  handleClick() {

  }
  render () {
    if (!this.props.expression) return null
    return (
      <div style={{marginTop: 50}}>
        <SectionHeader text="Expression"/>
        {
          _.map(this.props.expression.genes, (gene) => {
            return (
              <button
                key={gene}
                onClick={this.handleClick.bind(this)}
                style={{marginRight: 10}}>
                {gene}
              </button>
            )
          })
        }
      </div>
    )
  }
};

export default Expression;
