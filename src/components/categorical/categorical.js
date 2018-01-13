import React from "react";
import _ from "lodash";
import { connect } from "react-redux";

import * as globals from "../../globals";
import styles from "./categorical.css";
import SectionHeader from "../framework/sectionHeader"
import Value from "./value";
import { alphabeticallySortedValues } from "./util";

@connect((state) => {
  return {
    colorAccessor: state.controls.colorAccessor 
  }
})
class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isChecked: true
    };
  }

  handleColorChange() {
    this.props.dispatch({
      type: "color by categorical metadata",
      colorAccessor: this.props.metadataField
    })
  }
  toggleAll() {
    this.props.dispatch({
      type: "categorical metadata filter all of these",
      metadataField: this.props.metadataField
    })
    this.setState({isChecked: true})
  }
  toggleNone() {
    this.props.dispatch({
      type: "categorical metadata filter none of these",
      metadataField: this.props.metadataField,
      value: this.props.value
    })
    this.setState({isChecked: false})
  }
  render() {

    return (
      <div style={{
          // display: "flex",
          // alignItems: "baseline",
          maxWidth: globals.maxControlsWidth,
        }}>
        <div style={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "baseline"
        }}>
          <p style={{
              // flexShrink: 0,
              fontWeight: 500,
              // textAlign: "right",
              // fontFamily: globals.accentFont,
              // fontStyle: "italic",
              marginRight: 10,
            }}>
            {this.props.metadataField}
            <input
            onChange={this.state.isChecked ? this.toggleNone.bind(this) : this.toggleAll.bind(this)}
            checked={this.state.isChecked}
            type="checkbox"/>
            <span
              onClick={this.handleColorChange.bind(this)}
              style={{
                fontSize: 16,
                marginLeft: 4,
                position: "relative",
                top: 2,
                cursor: "pointer",
              }}>
              🖌️
            </span>
          </p>
        </div>
        <div>
          {
            _.map(alphabeticallySortedValues(this.props.values), (v, i) => {
              return (
                <Value
                  key={v}
                  metadataField={this.props.metadataField}
                  count={this.props.values[v]}
                  value={v}
                  i={i} />
              )
            })
          }
        </div>
      </div>
    )
  }
}

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

  render () {
    if (!this.props.ranges) return null

    return (
      <div style={{
        width: 310,
        marginRight: 40,
        paddingRight: 20,
        flexShrink: 0,
        height: 700,
        overflow: "auto",
      }}>

        {
          _.map(this.props.ranges, (value, key) => {
            const isColorField = key.includes("color") || key.includes("Color");
            if (
              value.options &&
              key !== "CellName" &&
              !isColorField
            ) {
              return (
                <Category
                  key={key}
                  metadataField={key}
                  values={value.options}/>
              )
            }
          })
        }
      </div>
    )
  }
};

export default Categories;

/*
<SectionHeader text="Categorical Metadata"/>
  [on off] toggle hide deselected filters (shows a menu vs shows what you ordered in compact/narrative form. fold out animation.)

  <p> Sort by: Alphabetical / Count || Show counts: true / false || Collapse: all / none</p>


  <p> <button> Field name [initial state] [create 2d graph]             [explain part of 2d graph]    </button> </p>
  <p> <button> Field name [initial state] [create hypothesis]           [validate hypothesis]         </button> </p>
  <p> <button> Field name [total]         [selected in current filters] [selected in graph selection] </button> </p>
  <p> <button> Tumor      [400]           [40]                          [4]                           </button> </p>

  might be interesting to show them as an 👁 icon, or with a slash through it, to allow for visible or hidden state

*/


/*

  Each category has a color associated with it - ie., color by location should show up on these buttons,
  Does colorby live here? Like a global, with the category labels at the top? If not, we have to
  duplicate the category names somewhere else - but then, not all (like Cluster_2d_color) won't be listed here.

*/
