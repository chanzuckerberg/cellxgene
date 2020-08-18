/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import { connect } from "react-redux";
import HistogramBrush from "../brushableHistogram";

@connect((state) => ({
  schema: state.annoMatrix?.schema,
}))
class Continuous extends React.PureComponent {
  render() {
    /* initial value for iterator to simulate index, ranges is an object */
    const { schema } = this.props;
    if (!schema) return null;
    const obsIndex = schema.annotations.obs.index;
    const allContinuousNames = schema.annotations.obs.columns
      .filter((col) => col.type === "int32" || col.type === "float32")
      .filter((col) => col.name !== obsIndex)
      .filter((col) => !col.writable) // skip user annotations - they will be treated as categorical
      .map((col) => col.name);

    return (
      <div>
        {allContinuousNames.map((key, zebra) => (
          <HistogramBrush key={key} field={key} isObs zebra={zebra % 2 === 0} />
        ))}
      </div>
    );
  }
}

export default Continuous;
