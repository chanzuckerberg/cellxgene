import React from "react";
import { connect } from "react-redux";
import { InputGroup } from "@blueprintjs/core";

const VanillaInput = props => {
  const { text, handleTextChange, inputProps } = props;
  return (
    <InputGroup
      {...inputProps}
      autoFocus
      value={text}
      intent="none"
      onChange={e => handleTextChange?.(e.target.value)}
      leftIcon="tag"
    />
  );
};

@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  world: state.world,
  ontology: state.ontology?.terms,
  ontologyLoading: state.ontology?.loading
}))
class AnnoInputs extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const { handleTextChange, text, ...restProps } = this.props;
    return (
      <div>
        <VanillaInput
          {...restProps}
          text={text}
          handleTextChange={handleTextChange}
        />
      </div>
    );
  }
}

export default AnnoInputs;
