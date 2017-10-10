import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";

import cells from "../../../data/GBM_metadata.js";
import createCategoryCounts from "./createCategoryCounts";

const Button = ({category, c, i}) => (
  <button
    key={i}
    style={{
      border: i > -1 ? "none" : `1px solid ${globals.lightgrey}`,
      background: i > -1 ? globals.hcaBlue : "none",
      color: i > -1 ? "white" : "black",
      padding: "7px 10px",
      marginRight: 20,
      borderRadius: 3,
      cursor: "pointer"
    }}>
    {c}
    <span
      style={{
        fontWeight: i > -1 ? globals.bolder : "auto",
        background: i > -1 ? "white" : "none",
        color: "black",
        marginLeft: 10,
        padding: "2px 4px",
        borderRadius: 3
      }}>
      {category[c]}
    </span>
  </button>
)

const Category = ({category, title}) => (
  <div>
    <p>{title}</p>
    { _.map(Object.keys(category), (c, i) => <Button key={i} category={category} c={c} i={i} />) }
  </div>
)

const Categorical = () => {

  return (
    <div style={{marginTop: 50}}>
      <h3>Categorical Metadata</h3>
      {
        _.map(createCategoryCounts(cells),
        (value, key) => {
          return <Category key={key} title={key} category={value} />
        })
      }
    </div>
  )
        <SectionHeader text="Categorical Metadata"/>
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
