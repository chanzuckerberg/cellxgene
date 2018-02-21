import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import styles from "./expression.css";
import SectionHeader from "../framework/sectionHeader"
import actions from "../../actions";
import ReactAutocomplete from "react-autocomplete";
import getContrast from "font-color-contrast"; // https://www.npmjs.com/package/font-color-contrast

class HeatmapSquare extends React.Component {
  constructor (props) {
    super(props)
    this.state = {
      value: '',
    }
  }
  render() {
    const contrastColor = getContrast(
      this.props.backgroundColor
        .substring(4, this.props.backgroundColor.length-1)
        .replace(/ /g, "")
        .split(",")
      )
    return (
      <p style={{
        padding: "12px 6px",
        textAlign: "center",
        color: contrastColor,
        width: 40,
        fontSize: 12,
        margin: 0,
        backgroundColor: this.props.backgroundColor,
      }}>
        {this.props.text}
      </p>
    )
  }
}

/**********************************
***********************************
***********************************
              Row
***********************************
***********************************
**********************************/
@connect((state) => {
  
  return {
    scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
    scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
    colorAccessor: state.controls.colorAccessor,
  }
})
class HeatmapRow extends React.Component {
  constructor (props) {
    super(props)
    this.state = {
      value: '',
    }
  }
  handleGeneColorScaleClick(gene) {
    return () => {
      this.props.dispatch(
        actions.requestSingleGeneExpressionCountsForColoringPOST(this.props.gene)
      )
    }
  }
  handleSetGeneAsScatterplotX(gene) {
    return () => {
      this.props.dispatch({
        type: "set scatterplot x",
        data: this.props.gene
      })
    }
  }
  handleSetGeneAsScatterplotY(gene) {
    return () => {
      this.props.dispatch({
        type: "set scatterplot y",
        data: this.props.gene
      })
    }
  }
  render() {
    return (
      <div style={{
        width: 220,
        display: "flex",
        justifyContent: "flex-start",
        alignItems: "baseline",
      }}>
        <div style={{width: 200}}>
        <span
          onClick={this.handleSetGeneAsScatterplotX(this.props.gene).bind(this)}
          style={{
            fontSize: 16,
            color: "#0074D9", // http://clrs.cc/
            cursor: "pointer",
            position: "relative",
            top: 1,
            fontWeight: 700,
            marginRight: 6,
            border: this.props.scatterplotXXaccessor === this.props.gene ? "1px solid black" : "none",
          }}>X</span>
          <span
            onClick={this.handleSetGeneAsScatterplotY(this.props.gene).bind(this)}
            style={{
              fontSize: 16,
              color: "#0074D9", // http://clrs.cc/
              cursor: "pointer",
              position: "relative",
              top: 1,
              fontWeight: 700,
              marginRight: 6,
              border: this.props.scatterplotYYaccessor === this.props.gene ? "1px solid black" : "none",
            }}>Y</span>
          <span
            onClick={this.handleGeneColorScaleClick(this.props.gene).bind(this)}
            style={{
              fontSize: 16,
              cursor: "pointer",
              position: "relative",
              top: 2,
              marginRight: 6,
              border: this.props.colorAccessor === this.props.gene ? "1px solid black" : "none",
            }}>🖌️</span>
          <span
            style={{
              fontSize: 14
            }}>
            {this.props.gene}
          </span>
        </div>
        <HeatmapSquare
          backgroundColor={this.props.greyColorScale(this.props.set1exp)}
          text={this.props.set1exp}/>
        <HeatmapSquare
          backgroundColor={this.props.greyColorScale(this.props.set2exp)}
          text={this.props.set2exp}/>
      </div>
    )
  }
}

/**********************************
***********************************
***********************************
        DiffExp Heatmap
***********************************
***********************************
**********************************/

@connect((state) => {
  return {
    differential: state.differential,
    allGeneNames: state.controls.allGeneNames
  }
})
class Heatmap extends React.Component {
  constructor (props) {
    super(props)
    this.state = {
      value: '',
    }
  }
  render() {
    if (!this.props.differential.diffExp) return null

    const topGenesForCellSet1 = this.props.differential.diffExp.data.celllist1;
    const topGenesForCellSet2 = this.props.differential.diffExp.data.celllist2;

    const extent = d3.extent(
      _.union(
        topGenesForCellSet1.mean_expression_cellset1,
        topGenesForCellSet1.mean_expression_cellset2,
        topGenesForCellSet2.mean_expression_cellset1,
        topGenesForCellSet2.mean_expression_cellset2
      )
    )

    const greyColorScale = d3.scaleSequential()
        .domain(extent)
        .interpolator(d3.interpolateGreys);

    return (
      <div>
        Color by any gene:
        <ReactAutocomplete
          items={this.props.allGeneNames}
          shouldItemRender={(item, value) => item.toLowerCase().indexOf(value.toLowerCase()) > -1}
          getItemValue={item => item}
          renderItem={(item, highlighted) =>
            <div
              key={item}
              style={{ backgroundColor: highlighted ? '#eee' : 'transparent'}}
            >
              {item}
            </div>
          }
          value={this.state.value}
          onChange={e => this.setState({ value: e.target.value })}
          onSelect={(value) => {
            this.setState({ value });
            this.props.dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(value))
          }}
        />
        <div style={{
          display: "flex",
          justifyContent: "space-between",
          width: 165,
          fontWeight: 700,
        }}>
          <p style={{width: 100}}>Gene</p>
          <p>1</p>
          <p>2</p>
        </div>
        {
          topGenesForCellSet1.topgenes.map((gene, i) => {
            return <HeatmapRow
              key={gene}
              gene={gene}
              greyColorScale={greyColorScale}
              set1exp={Math.floor(topGenesForCellSet1.mean_expression_cellset1[i])}
              set2exp={Math.floor(topGenesForCellSet1.mean_expression_cellset2[i])}
              />
          })
        }
        {
          topGenesForCellSet2.topgenes.map((gene, i) => {
            return <HeatmapRow
              key={gene}
              gene={gene}
              greyColorScale={greyColorScale}
              set1exp={Math.floor(topGenesForCellSet2.mean_expression_cellset1[i])}
              set2exp={Math.floor(topGenesForCellSet2.mean_expression_cellset2[i])}
              />
          })
        }
      </div>
    )
  }
}

export default Heatmap;
