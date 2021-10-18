import React from "react";
import { connect } from "react-redux";
import {
  ControlGroup,
} from "@blueprintjs/core";
import ParameterInput from "./parameterinput";
import { ControlsHelpers } from "../../util/stateManager";
import DefaultsButton from "./defaultsio";

@connect((state) => ({
  reembedParams: state.reembedParameters,
  annoMatrix: state.annoMatrix
}))
class BatchPanel extends React.PureComponent {
  
  render() {
    const { reembedParams, annoMatrix, dispatch } = this.props;
    const allCategoryNames = ControlsHelpers.selectableCategoryNames(
      annoMatrix.schema
    ).sort();
    const disabled = !reembedParams.doBatch
    let panel;
    switch (reembedParams.batchMethod) { 
      case "BBKNN": {
        panel = (
          <div style={{"paddingLeft":"10px"}}>
            <ParameterInput 
              min={1}
              label="neighbors_within_batch"
              param="bbknnNeighborsWithinBatch"
            />                   
          </div>
        );
        break;
      } case "Scanorama": {
        panel = (
          <div style={{"paddingLeft":"10px"}}>
            <ControlGroup fill={true} vertical={false}>
              <ParameterInput 
                min={1}
                label="knn"
                param="scanoramaKnn"
              />     
              <ParameterInput 
                min={0}
                label="sigma"
                param="scanoramaSigma"
              />      
            </ControlGroup>  
            <ControlGroup fill={true} vertical={false}>             
              <ParameterInput 
                min={0}
                max={1}
                label="alpha"
                param="scanoramaAlpha"
              />  
              <ParameterInput 
                min={0}
                label="batch_size"
                param="scanoramaBatchSize"
              />                                                   
            </ControlGroup>
          </div>
        );
        break;
      } case "Harmony": {
        panel = (
          <div style={{"paddingLeft":"10px"}}>
            <ControlGroup fill={true} vertical={false}>

            </ControlGroup>
          </div>
        );
        break;
      } default: {
        panel = null;
      }
    }
    panel = reembedParams.doBatch ? panel : null;
    return (
      <div>
        <ControlGroup fill={true} vertical={false}>
          <ParameterInput 
            label="Batch correct?"
            param="doBatch"
          />
          <ParameterInput 
            disabled={disabled}
            label="Method"
            param="batchMethod"
            options={["BBKNN","Harmony","Scanorama"]}
          />      
          <ParameterInput 
            disabled={disabled}        
            label="Batch key"
            param="batchKey"
            options={allCategoryNames}
          />    
        </ControlGroup>  
      {panel}      
    </div>
    );
  }
}

export default BatchPanel;
