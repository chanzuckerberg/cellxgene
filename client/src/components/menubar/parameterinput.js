import React, { useEffect } from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  NumericInput,
  Label,
  Checkbox,
  MenuItem
} from "@blueprintjs/core";
import { defaultPrepParams } from "../../reducers/reembed";
import { Select } from "@blueprintjs/select";
import { isNumber } from "lodash";

@connect((state) => ({
  reembedParams: state.reembedParameters,
}))
class ParameterInput extends React.PureComponent {
  constructor(props) {
    super(props);
    const { reembedParams, param } = props;
    let params = reembedParams;
    if (param in defaultPrepParams) {
      if (reembedParams.batchPrepKey !== "" && reembedParams.batchPrepLabel !== "") {
        params = reembedParams.batchPrepParams[reembedParams.batchPrepKey][reembedParams.batchPrepLabel]
      }      
    }    
    this.state = {
      value: params[param].toString(),
      refresher: false
    }
  }
  componentDidUpdate(prevProps) {
    const { reembedParams, param } = this.props;
    if (prevProps.reembedParams.batchPrepLabel !== reembedParams.batchPrepLabel){
      this.setState({
        refresher: true
      })
    }else if (reembedParams[param] !== prevProps.reembedParams[param] && isNumber(reembedParams[param])) {
      this.setState({ 
        value: reembedParams[param].toString()
      });      
    }
  }

  clamp = (num, min=Number.POSITIVE_INFINITY, max=Number.NEGATIVE_INFINITY) => {
    return Math.min(Math.max(num, min), max);
  };      
  render() {
    const { dispatch, param, label, min, max, reembedParams } = this.props;    
    let { value, refresher } = this.state;
    let params = reembedParams;

    if (param in defaultPrepParams) {
      if (reembedParams.batchPrepKey !== "" && reembedParams.batchPrepLabel !== "") {
        params = reembedParams.batchPrepParams[reembedParams.batchPrepKey][reembedParams.batchPrepLabel]
      }      
    }
    switch (typeof params[param]) {
      case "boolean": {
        return (
          <div>
            <Checkbox checked={params[param]} label={label} style={{"paddingTop":"10px"}}
              onChange={() => {
                dispatch({
                  type: "reembed: set parameter",
                  key: param,
                  value: !params[param]
                  })
                }
              } 
            /> 
          </div>          
        )
      } case "string": {
        const { disabled, options } = this.props;
        return (
          <div style={{"paddingTop":"5px"}}>
            <Select
            disabled={disabled}
            items={
              options
            }
            filterable={false}
            itemRenderer={(d, { handleClick }) => {
              return (
                <MenuItem
                  onClick={handleClick}
                  key={d}
                  text={d}
                />
              );
            }}
            onItemSelect={(d) => {
              dispatch({
                type: "reembed: set parameter",
                key: param,
                value: d
                }) 
            }}
          >
            <AnchorButton
              disabled={disabled}
              text={`${label}: ${params[param]}`}
              rightIcon="double-caret-vertical"
            />
          </Select>
        </div>
        );
      } default: {
        const { disabled } = this.props;
        return (
          <Label>
            {label}
            <NumericInput
              disabled={disabled}
              allowNumericCharactersOnly={true}
              placeholder={label}
              value={refresher ? params[param] : value}
              min={min}
              max={max}
              minorStepSize={0.001}
              onValueChange={
                (_valueAsNumber, valueAsString) => {
                  let val = (valueAsString.charAt(0)==="0" && valueAsString.charAt(1)!=="." 
                    ? valueAsString.substr(1) : valueAsString
                  );
                  val =  val!=="" && parseFloat(val)>=min ? val : min.toString();
                  const clamped = parseFloat(val) < min || parseFloat(val) > max;
                  val = clamped ? this.clamp(parseFloat(val),min,max).toString() : val;
                  this.setState({value: val, refresher: false})
                  dispatch({
                    type: "reembed: set parameter",
                    key: param,
                    value: parseFloat(val)
                  })
                }
              }
            />
          </Label> 
        );
      }
    }    
  }
}
export default ParameterInput;
