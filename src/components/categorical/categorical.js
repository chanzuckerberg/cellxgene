import React from "react";
import _ from "lodash";
import { connect } from "react-redux";

import * as globals from "../../globals";
import styles from "./categorical.css";
import SectionHeader from "../framework/sectionHeader"
import createCategoryCounts from "./createCategoryCounts";
import Value from "./value";
import { alphabeticallySortedValues } from "./util";

const Category = ({metadataField, values, toggleOn, toggleOff}) => (
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
            <Value
              key={v}
              toggleOn={toggleOn}
              toggleOff={toggleOff}
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
  toggleOn(metadataField, value) {
    return () => {
      this.props.dispatch({
        type: "categorical metadata filter selected",
        metadataField,
        value
      })
    }
  }
  toggleOff(metadataField, value) {
    return () => {
      this.props.dispatch({
        type: "categorical metadata filter deselected",
        metadataField,
        value
      })
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
                  toggleOn={this.toggleOn.bind(this)}
                  toggleOff={this.toggleOff.bind(this)}
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
