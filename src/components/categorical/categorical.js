import React from "react";
import _ from "lodash";
import { connect } from "react-redux";

import * as globals from "../../globals";
import styles from "./categorical.css";
import SectionHeader from "../framework/sectionHeader"
import createCategoryCounts from "./createCategoryCounts";
import { alphabeticallySortedValues } from "../../util";

import cells from "../../../data/GBM_metadata.js";

class Button extends React.Component {
  render () {
    return (
      <button
        key={this.props.i}
        onClick={this.props.handleClick(this.props.category, this.props.value)}
        style={{
          border: "none",
          background: "none",
          color: "black",
          padding: "7px 10px",
          marginRight: 3,
          cursor: "pointer"
        }}>
        {this.props.value}
        <span
          style={{
            fontFamily: globals.accentFont,
            fontStyle: "italic",
            color: "black",
            marginLeft: 4,
            borderRadius: 3
          }}>
          {this.props.count}
        </span>
      </button>
    )
  }
}


const Category = ({category, values, handleClick}) => (
  <div style={{
      display: "flex",
      alignItems: "baseline",
      maxWidth: globals.maxParagraphWidth,
    }}>
    <p style={{
        flexShrink: 0,
        width: 100,
        textAlign: "right",
        marginRight: 20
      }}>{category}:</p>
    <div>
      {
        _.map(alphabeticallySortedValues(values), (v, i) => {
          return (
            <Button
              key={v}
              handleClick={handleClick}
              category={category}
              count={values[v]}
              value={v}
              i={i} />
          )
        })
      }
    </div>
  </div>
)

@connect((state) => {
  return {
    foo123: state
  }
})
class Categorical extends React.Component {
  constructor(props) {
    super(props);
    this.state = {

    };
  }
  handleClick(category, value) {
    return () => {
      this.props.dispatch({type: "category changed", data: {category, value}})
    }
  }
  render () {
    return (
      <div style={{marginTop: 50}}>
        <SectionHeader text="Categorical Metadata"/>
        <button onClick={this.handleClick("whee", "wooo")}></button>
        <p> Sort by: Alphabetical / Count || Show counts: true / false || Collapse: all / none</p>
        {
          _.map(createCategoryCounts(cells),
          (value, key) => {
            return (
              <Category
                handleClick={this.handleClick.bind(this)}
                key={key}
                category={key}
                values={value} />
            )
          })
        }
      </div>
    )
  }
};

export default Categorical;

/*

  [on off] toggle hide deselected filters (shows a menu vs shows what you ordered in compact/narrative form. fold out animation.)

  <p> <button> Field name [initial state] [create 2d graph]             [explain part of 2d graph]    </button> </p>
  <p> <button> Field name [initial state] [create hypothesis]           [validate hypothesis]         </button> </p>
  <p> <button> Field name [total]         [selected in current filters] [selected in graph selection] </button> </p>
  <p> <button> Tumor      [400]           [40]                          [4]                           </button> </p>

*/


/*

  Each category has a color associated with it - ie., color by location should show up on these buttons,
  Does colorby live here? Like a global, with the category labels at the top? If not, we have to
  duplicate the category names somewhere else - but then, not all (like Cluster_2d_color) won't be listed here.

*/
