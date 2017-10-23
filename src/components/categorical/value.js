
import { connect } from "react-redux";
import React from "react";
import * as globals from "../../globals";

@connect((state) => {
  return {
    selectedMetadata: state.selectedMetadata
  }
})
class CategoryValue extends React.Component {
  render () {

    /* has this value for this category already been selected? */
    let selected = false;
    let fontWeight;

    if (!this.props.selectedMetadata) { /* there is no selected metadata */
      selected = false;
    } else if (
      this.props.selectedMetadata[this.props.metadataField] && /* the key {Location: []} is present  */
      this.props.selectedMetadata[this.props.metadataField].indexOf(this.props.value) > -1 /* "Tumor" exists in {Location: ["Tumor"]}  */
    ) {
      selected = true;
    }

    return (
      <div
        key={this.props.i}
        onClick={selected ? this.props.toggleOff(this.props.metadataField, this.props.value) : this.props.toggleOn(this.props.metadataField, this.props.value)}
        style={{
          cursor: "pointer",
          display: "flex",
          fontWeight: selected ? 700 : 400,
        }}>
        <p style={{
            width: 200,
            flexShrink: 0,
            margin: 0,
            lineHeight: "1em"
          }}>
          {this.props.value}
        </p>
        <p style={{
            margin: 0,
            lineHeight: "1em"
          }}>
          {this.props.count}
        </p>
      </div>
    )
  }
}

export default CategoryValue;
