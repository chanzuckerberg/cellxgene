import React from "react";
import _ from "lodash";
import { connect } from "react-redux";

import * as globals from "../../globals";
import styles from "./categorical.css";
import SectionHeader from "../framework/sectionHeader"
import createCategoryCounts from "./createCategoryCounts";
import { alphabeticallySortedValues } from "./util";

@connect((state) => {
  return {
    selectedMetadata: state.selectedMetadata
  }
})
class Button extends React.Component {
  render () {

    /* has this value for this category already been selected? */
    let fontWeight;
    if (!this.props.selectedMetadata) { /* there is no selected metadata */
      fontWeight = 400;
    } else if (
      this.props.selectedMetadata[this.props.metadataField] && /* the key {Location: []} is present  */
      this.props.selectedMetadata[this.props.metadataField].indexOf(this.props.value) > -1 /* "Tumor" exists in {Location: ["Tumor"]}  */
    ) {
      fontWeight = 700
    }

    return (
      <div
        key={this.props.i}
        onClick={this.props.handleClick(this.props.metadataField, this.props.value)}
        style={{
          cursor: "pointer",
          display: "flex",
          fontWeight,
        }}>
        <div style={{
            width: 200,
            flexShrink: 0,
          }}>
          {this.props.value}
        </div>
        <span>
          {this.props.count}
        </span>
      </div>
    )
  }
}

const Category = ({metadataField, values, handleClick}) => (
  <div style={{
      display: "flex",
      alignItems: "baseline",
      maxWidth: globals.maxControlsWidth,
    }}>
    <p style={{
        flexShrink: 0,
        width: 100,
        textAlign: "right",
        fontFamily: globals.accentFont,
        fontStyle: "italic",
        marginRight: 20,
      }}>{metadataField}:</p>
    <div>
      {
        _.map(alphabeticallySortedValues(values), (v, i) => {
          return (
            <Button
              key={v}
              handleClick={handleClick}
              metadataField={metadataField}
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

  const ranges = state.cells.cells && state.cells.cells.data.ranges ? state.cells.cells.data.ranges : null;

  return {
    ranges
  }
})
class Categories extends React.Component {
  constructor(props) {
    super(props);
    this.state = {

    };
  }
  handleClick(metadataField, value) {
    return () => {
      this.props.dispatch({type: "categorical metadata filter selected", metadataField, value})
    }
  }
  render () {
    if (!this.props.ranges) return null
    return (
      <div style={{marginTop: 50}}>
        <SectionHeader text="Categorical Metadata"/>
        {
          _.map(this.props.ranges, (value, key) => {
            if (
              value.options &&
              key !== "CellName"
            ) {
              return (
                <Category
                  handleClick={this.handleClick.bind(this)}
                  key={key}
                  metadataField={key}
                  values={value.options}/>
              )
            } else {
              return null
            }
          })
        }
      </div>
    )
  }
};

export default Categories;

/*

  [on off] toggle hide deselected filters (shows a menu vs shows what you ordered in compact/narrative form. fold out animation.)

  <p> Sort by: Alphabetical / Count || Show counts: true / false || Collapse: all / none</p>


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
