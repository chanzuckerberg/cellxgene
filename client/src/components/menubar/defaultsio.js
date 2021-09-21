import React from "react";
import {
  AnchorButton,
  ControlGroup,
} from "@blueprintjs/core";
import { reembedParamsFetch } from "../../actions";
import actions from "../../actions";

class DefaultsButton extends React.PureComponent {
  
  render() {
    const { dispatch } = this.props;
    return (     
    <ControlGroup fill={true} vertical={false}>
      <AnchorButton
        onClick={async () => {
          await reembedParamsFetch(dispatch)
          }
        } 
        text={`Load defaults`}
        small outlined
      />           
      <AnchorButton
        onClick={async () => {
          dispatch(actions.saveReembedParametersAction())
        }}
        text={`Save as default`}
        small
        outlined
      />                         
    </ControlGroup> 
    );
  }
}

export default DefaultsButton;
