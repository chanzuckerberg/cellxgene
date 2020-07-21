import React from "react";
import { connect } from "react-redux";

@connect((state) => ({
  layoutChoice: state.layoutChoice,
}))
class Embedding extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const { layoutChoice } = this.props;
    return (
      <div
        id="embedding"
        data-testclass="embedding"
        style={{
          position: "absolute",
          display: "inherit",
          left: 8,
          bottom: 8,
          zIndex: 1,
          fontWeight: 700,
        }}
      >
        {layoutChoice?.current}
      </div>
    );
  }
}

export default Embedding;
