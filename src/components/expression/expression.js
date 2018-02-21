import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import styles from "./expression.css";
import actions from "../../actions";
import Heatmap from "./diffExpHeatmap";
// <p> Top genes expressed by cellset 2 </p>
// <p> Gene {"............"} Exp in set 1 {"............"} Exp in set 2 </p>

@connect((state) => {
  return {
    currentCellSelection: state.controls.currentCellSelection,
    differential: state.differential
  }
})
class Expression extends React.Component {
  constructor(props) {
    super(props);
    this.state = {

    };
  }
  handleClick(gene) {
    return () => {
      this.props.dispatch({
        type: "color by expression",
        gene: gene,
      });
    }
  }
  set1() {
    const set = [];

    _.each(this.props.currentCellSelection, (cell) => {
      if (cell["__selected__"]) {
        set.push(cell.CellName)
      }
    })

    this.props.dispatch({
      type: "store current cell selection as differential set 1",
      data: set
    })
  }
  set2() {
    const set = [];

    _.each(this.props.currentCellSelection, (cell) => {
      if (cell["__selected__"]) {
        set.push(cell.CellName)
      }
    })

    this.props.dispatch({
      type: "store current cell selection as differential set 2",
      data: set
    })
  }
  computeDiffExp() {
    this.props.dispatch(
      actions.requestDifferentialExpression(
        this.props.differential.celllist1,
        this.props.differential.celllist2
      )
    )
  }
  render () {
    if (
      !this.props.differential
    ) {
      return null
    }
    return (
      <div>
        <div style={{
          marginBottom: 20,
          width: 500,
        }}>
          <button
            style={{
              width: 200,
              fontSize: 14,
              margin: 5,
              fontWeight: 400,
              color: "#5A5F63",
              padding: "10px 20px",
              backgroundColor: "#C5D1D8",
              border: "none",
              borderRadius: 3,
              cursor: "pointer",

            }}
            onClick={this.set1.bind(this)}>
            {
              this.props.differential.celllist1 ?
              this.props.differential.celllist1.length + " cells selected" :
              "Store current cell selection as 'set 1'"
            }
          </button>
          <button
            style={{
              width: 200,
              fontSize: 14,
              margin: 5,
              fontWeight: 400,
              color: "#5A5F63",
              padding: "10px 20px",
              backgroundColor: "#C5D1D8",
              border: "none",
              borderRadius: 3,
              cursor: "pointer",
            }}
            onClick={this.set2.bind(this)}>
            {
              this.props.differential.celllist2 ?
              this.props.differential.celllist2.length + " cells selected" :
              "Store current cell selection as 'set 2'"
            }
          </button>
          </div>
        <div>
          <button
            style={{
              fontSize: 14,
              margin: 5,
              fontWeight: 400,
              color: "#41633C",
              padding: "10px 20px",
              backgroundColor: "#B9D1B5",
              border: "none",
              borderRadius: 3,
              cursor: "pointer",
            }}
            onClick={this.computeDiffExp.bind(this)}>
            Compute differential expression
          </button>
        </div>
        <Heatmap/>
      </div>
    )
  }
};

export default Expression;
