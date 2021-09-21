import React from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  Collapse,
  ControlGroup,
} from "@blueprintjs/core";
import ParameterInput from "./parameterinput";
import DefaultsButton from "./defaultsio";

@connect((state) => ({
  reembedParams: state.reembedParameters,
  annoMatrix: state.annoMatrix
}))
class PrepPanel extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      cfshown: false,
      gfshown: false,
      hvgshown: false,
    };
  }  
  componentDidUpdate(prevProps) {
    const { reembedParams } = this.props;
    if (reembedParams.doPreprocess !== prevProps.reembedParams.doPreprocess && !reembedParams.doPreprocess) {
      this.setState({ 
        hvgshown: false,
        cfshown: false,
        gfshown: false,        
      });
    }
  }

  render() {
    const {
      cfshown, gfshown, hvgshown
    } = this.state;
    const { reembedParams, annoMatrix, dispatch} = this.props;
    
    const disabled = !reembedParams.doPreprocess
    return (
      <div>
      <DefaultsButton dispatch={dispatch}/>
      <ControlGroup fill={true} vertical={false}>
        <ParameterInput 
          label="Preprocess?"
          param="doPreprocess"
        />
        <ParameterInput 
          label="SAM?"
          param="doSAM"
        />                   
      </ControlGroup>    

      <AnchorButton
        onClick={() => {
          this.setState({ 
            hvgshown: false,
            cfshown: !this.state.cfshown,
            gfshown: false,
          });
        }}
        text={`Cell filtering`}
        fill outlined
        rightIcon={cfshown ? "chevron-down" : "chevron-right"} small
        disabled = {true}
      />                    
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={cfshown}>
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
          });
        }}
        text={`Gene filtering`}
        fill outlined
        rightIcon={gfshown ? "chevron-down" : "chevron-right"} small
        disabled = {disabled}
      />   
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={gfshown}>
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
            hvgshown: !this.state.hvgshown,
            cfshown: false,
            gfshown: false,
          });
        }}
        text={`Highly variable gene selection`}
        fill outlined
        rightIcon={hvgshown ? "chevron-down" : "chevron-right"} small
        disabled = {disabled || reembedParams.doSAM}
      />   
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={hvgshown}>   
          <ControlGroup fill={true} vertical={false}>
            <ParameterInput
              min={0}
              disabled={reembedParams.doSAM}
              max={annoMatrix.nVar}
              label="n_top_genes"
              param="nTopGenesHVG"
            />  
            <ParameterInput
              min={1}
              disabled={reembedParams.doSAM}
              label="n_bins"
              param="nBinsHVG"
            />        
          </ControlGroup>                    
        </Collapse>  
      </div>                  
    </div>
    );
  }
}

export default PrepPanel;
