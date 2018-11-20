// jshint esversion: 6
import React from "react";
import _ from "lodash";

@connect()
class Occupancy extends React.Component {
  render() {
    const {
      categoricalAsBooleansMap,
      metadataField,
      count,
      value,
      colorAccessor,
      colorScale,
      schema
    } = this.props;

    return (
      <svg
        style={{
          marginLeft: 5,
          width: 11,
          height: 11,
          backgroundColor: "red"
        }}
      />
    );
  }
}

export default Occupancy;
