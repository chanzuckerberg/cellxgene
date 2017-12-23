import { connect } from "react-redux";
import React from "react";
import * as globals from "../../globals";
import actions from "../../actions";

@connect((state) => {
  return {
    selectedMetadata: state.selectedMetadata
  }
})
class CategoryValue extends React.Component {

  handleChange() {
    console.log("dispatch", this.props.metadataField, this.props.value)
  }

  render () {

    /* has this value for this category already been selected? */
    let selected = false;

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
        style={{
          display: "flex",
          fontWeight: selected ? 700 : 400,
        }}>
        <p style={{
            width: 200,
            flexShrink: 0,
            margin: 0,
            lineHeight: "1em"
          }}>
          <input onChange={this.handleChange.bind(this)} checked={true} type="checkbox"/> {this.props.value}
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

// onClick={selected ? this.toggleOff.bind(this) : this.toggleOn.bind(this)}
// toggleOn() {
//   this.props.dispatch(
//     actions.attemptCategoricalMetadataSelection(
//       this.props.metadataField,
//       this.props.value
//     ))
// }
// toggleOff() {
//   this.props.dispatch(
//     actions.attemptCategoricalMetadataDeselection(
//       this.props.metadataField,
//       this.props.value
//     ))
// }
