import React from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  Collapse,
  ControlGroup,
} from "@blueprintjs/core";
import ParameterInput from "./parameterinput";
import DefaultsButton from "./defaultsio";
import { ControlsHelpers } from "../../util/stateManager";

@connect((state) => ({
  reembedParams: state.reembedParameters,
  annoMatrix: state.annoMatrix,
}))
class PrepPanel extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      cfshown: false,
      gfshown: false,
      hvgshown: false,
      trshown: false
    };
  }  
  componentDidUpdate(prevProps) {
    const { reembedParams } = this.props;
    if (reembedParams.doPreprocess !== prevProps.reembedParams.doPreprocess && !reembedParams.doPreprocess) {
      this.setState({ 
        hvgshown: false,
        cfshown: false,
        gfshown: false,   
        trshown: false     
      });
    }
  }

  render() {
    const {
      cfshown, gfshown, hvgshown, trshown
    } = this.state;
    const { reembedParams, annoMatrix, dispatch} = this.props;
    const allCategoryNames = ControlsHelpers.selectableCategoryNames(
      annoMatrix.schema
    ).sort();
    let disabled = !reembedParams.doPreprocess
    const batchDisabled = !reembedParams.doBatchPrep
    let allBatchPrepLabels;
    let disabledBatchLabel = true;
    if (reembedParams.batchPrepKey !== ""){
      if (reembedParams.batchPrepLabel !== ""){
        disabled = !reembedParams.batchPrepParams[reembedParams.batchPrepKey][reembedParams.batchPrepLabel].doPreprocess
      }
      allBatchPrepLabels = annoMatrix.schema.annotations.obsByName?.[reembedParams.batchPrepKey]?.categories
      if(allBatchPrepLabels){
        allBatchPrepLabels = allBatchPrepLabels.filter(item => item !== "unassigned")
        disabledBatchLabel = false;
      }
    }
    allBatchPrepLabels = allBatchPrepLabels ?? [""]
    return (
      <div>
      <DefaultsButton dispatch={dispatch}/>       
      <ControlGroup fill={true} vertical={false}>     
        <ParameterInput 
          label="Batch preprocess?"
          param="doBatchPrep"
        />   
        <ParameterInput 
          disabled={batchDisabled}        
          label="Batch key"
          param="batchPrepKey"
          options={allCategoryNames}
        />
        <ParameterInput 
          disabled={batchDisabled || disabledBatchLabel}        
          label="Batch label"
          param="batchPrepLabel"
          options={allBatchPrepLabels}
        />              
      </ControlGroup>        
      <ControlGroup fill={true} vertical={false}>
        <ParameterInput 
          label="Preprocess?"
          param="doPreprocess"
        />           
        <ParameterInput
          label="Data layer"
          param="dataLayer"
          options={annoMatrix.schema.layers}
        />                   
      </ControlGroup>  
      <AnchorButton
        onClick={() => {
          this.setState({ 
            hvgshown: false,
            cfshown: !this.state.cfshown,
            gfshown: false,
            trshown: false
          });
        }}
        text={`Cell filtering`}
        fill outlined
        rightIcon={cfshown ? "chevron-down" : "chevron-right"} small
        disabled = {disabled}
      />                    
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={cfshown && !disabled}>
          <ControlGroup fill={true} vertical={false}>
            <ParameterInput
            min={0}
            label="min_counts"
            param="minCountsCF"
            />
            <ParameterInput
            min={0}
            label="min_genes"
            param="minGenesCF"
            />         
          </ControlGroup>
        </Collapse>
      </div>     
      <AnchorButton
        onClick={() => {
          this.setState({ 
            hvgshown: false,
            gfshown: !this.state.gfshown,
            cfshown: false,
            trshown: false
          });
        }}
        text={`Gene filtering`}
        fill outlined
        rightIcon={gfshown ? "chevron-down" : "chevron-right"} small
        disabled = {disabled}
      />   
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={gfshown && !disabled}>
          <ControlGroup fill={true} vertical={false}>
            <ParameterInput
              min={0}
              label="min_counts"
              param="minCountsGF"
              />        
            </ControlGroup>
            <ControlGroup fill={true} vertical={false}>
              <ParameterInput
              min={0}
              max={100}
              label="min_cells (%)"
              param="minCellsGF"
              />
              <ParameterInput
              min={0}
              max={100}
              label="max_cells (%)"
              param="maxCellsGF"
              />         
            </ControlGroup>          
        </Collapse>
      </div>
      <AnchorButton
        onClick={() => {
          this.setState({ 
            hvgshown: false,
            cfshown: false,
            gfshown: false,
            trshown: !this.state.trshown
          });
        }}
        text={`Transformation`}
        fill outlined
        rightIcon={trshown ? "chevron-down" : "chevron-right"} small
        disabled = {disabled}
      />   
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={trshown}>     
          <ControlGroup fill={true} vertical={false}>
          <ParameterInput
              label="Sum normalize?"
              param="sumNormalizeCells"
            />            
            <ParameterInput
              label="Log transform?"
              param="logTransform"
            />
          </ControlGroup> 
        </Collapse>       
      </div>                               
    </div>
    );
  }
}

export default PrepPanel;
