import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import styles from "./expression.css";
import SectionHeader from "../framework/sectionHeader"
import actions from "../../actions";

@connect((state) => {
  return {
    differential: state.differential
  }
})
class DiffExp extends React.Component {
  handleGeneColorScaleClick(gene) {
    return () => {
      this.props.dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(gene))
    }
  }
  render() {
    if (!this.props.differential.diffExp) return null

    const topGenesForCellSet1 = this.props.differential.diffExp.data.celllist1;
    const topGenesForCellSet2 = this.props.differential.diffExp.data.celllist2;
    const rectSide = 10;
    const rightMargin = 0;

    const extent = d3.extent(
      _.union(
        topGenesForCellSet1.mean_expression_cellset1,
        topGenesForCellSet1.mean_expression_cellset2,
        topGenesForCellSet2.mean_expression_cellset1,
        topGenesForCellSet2.mean_expression_cellset2
      )
    )

    console.log('diffexp', this.props)
    console.log('extent', extent)

    const greyColorScale = d3.scaleSequential()
        .domain(extent)
        .interpolator(d3.interpolateGreys);

    return (
      <div>
      <div style={{
        display: "flex",
        justifyContent: "space-between",
        width: 400,
        fontWeight: 700,
      }}>
        <p style={{width: 200}}>Gene</p>
        <p>Set 1 exp</p>
        <p>Set 2 exp</p>
      </div>
      {
        topGenesForCellSet1.topgenes.map((gene, i) => {
          return (
            <div key={gene} style={{
              width: 400,
              display: "flex",
              justifyContent: "space-between",
              alignItems: "baseline",
            }}>
              <div style={{width: 200}}>
                <span
                  onClick={this.handleGeneColorScaleClick(gene).bind(this)}
                  style={{
                    fontSize: 24,
                    cursor: "pointer",
                    position: "relative",
                    top: 5,
                    marginRight: 10
                  }}>🎨</span>
                <span>{gene}</span>
              </div>
              <p style={{
                padding: 20,
                color: extent[1] === topGenesForCellSet1.mean_expression_cellset1[i] ? "white" : "black",
                width: 60,
                margin: 0,
                backgroundColor: greyColorScale(Math.floor(topGenesForCellSet1.mean_expression_cellset1[i])),
              }}>
                {Math.floor(topGenesForCellSet1.mean_expression_cellset1[i])}
              </p>
              <p style={{
                padding: 20,
                color: extent[1] === topGenesForCellSet1.mean_expression_cellset2[i] ? "white" : "black",
                width: 60,
                margin: 0,
                backgroundColor: greyColorScale(Math.floor(topGenesForCellSet1.mean_expression_cellset2[i])),
              }}>
                {Math.floor(topGenesForCellSet1.mean_expression_cellset2[i])}
              </p>
            </div>
          )
        })
      }
      {
        topGenesForCellSet2.topgenes.map((gene, i) => {
          return (
            <div key={gene} style={{
              width: 400,
              display: "flex",
              justifyContent: "space-between",
              alignItems: "baseline",
            }}>
              <div style={{width: 200}}>
                <span
                  onClick={this.handleGeneColorScaleClick(gene).bind(this)}
                  style={{
                    fontSize: 24,
                    cursor: "pointer",
                    position: "relative",
                    top: 5,
                    marginRight: 10
                  }}>🎨</span>
                <span>{gene}</span>
              </div>
              <p style={{
                padding: 20,
                color: extent[1] === topGenesForCellSet2.mean_expression_cellset1[i] ? "white" : "black",
                width: 60,
                margin: 0,
                backgroundColor: greyColorScale(Math.floor(topGenesForCellSet2.mean_expression_cellset1[i])),
              }}>
                {Math.floor(topGenesForCellSet2.mean_expression_cellset1[i])}
              </p>
              <p style={{
                padding: 20,
                color: extent[1] === topGenesForCellSet2.mean_expression_cellset2[i] ? "white" : "black",
                width: 60,
                margin: 0,
                backgroundColor: greyColorScale(Math.floor(topGenesForCellSet2.mean_expression_cellset2[i])),
              }}>
                {Math.floor(topGenesForCellSet2.mean_expression_cellset2[i])}
              </p>
            </div>
          )
        })
      }
      </div>
    )
  }
}

// <p> Top genes expressed by cellset 2 </p>
// <p> Gene {"............"} Exp in set 1 {"............"} Exp in set 2 </p>


@connect((state) => {
  return {
    expression: state.expression.data,
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
      !this.props.expression ||
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
        <DiffExp/>
      </div>
    )
  }
};

export default Expression;


// <SectionHeader text="Differential Expression"/>
// <SectionHeader text="Color by expression"/>
// {
//   _.map(this.props.expression.genes, (gene) => {
//     return (
//       <button
//         key={gene}
//         onClick={this.handleClick(gene).bind(this)}
//         style={{marginRight: 10}}>
//         {gene}
//       </button>
//     )
//   })
// }


// Gene count:
// <input placeholder={"Default is 5"}/>

// <svg>
//   <g transform="translate(200, 40)">
//     {
//       _.map(topGenesForCellSet1.topgenes, (gene, i) => {
//         return (
//           <g key={gene}>
//             <rect
//               width={rectSide}
//               height={rectSide}
//               x={i * (rectSide + rightMargin)}
//               fill={"rgb(255,0,0)"}
//               />
//           <g>
//         )
//       })
//     }
//   </g>
// </svg>
