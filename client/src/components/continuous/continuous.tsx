/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import { connect } from "react-redux";
import HistogramBrush from "../brushableHistogram";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  schema: (state as any).annoMatrix?.schema,
}))
class Continuous extends React.PureComponent {
  render() {
    /* initial value for iterator to simulate index, ranges is an object */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
    const { schema } = this.props;
    if (!schema) return null;
    const obsIndex = schema.annotations.obs.index;
    const allContinuousNames = schema.annotations.obs.columns
      .filter((col: any) => col.type === "int32" || col.type === "float32")
      .filter((col: any) => col.name !== obsIndex)
      .filter((col: any) => !col.writable) // skip user annotations - they will be treated as categorical
      .map((col: any) => col.name);
    return (
      <div>
        {allContinuousNames.map((key: any, zebra: any) => (
          // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
          <HistogramBrush key={key} field={key} isObs zebra={zebra % 2 === 0} />
        ))}
      </div>
    );
  }
}

export default Continuous;
