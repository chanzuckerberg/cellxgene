import React from "react";
import { connect } from "react-redux";
import {
  AnchorButton,
  Collapse,
  ControlGroup,
  InputGroup
} from "@blueprintjs/core";
import ParameterInput from "./parameterinput";
import DefaultsButton from "./defaultsio";

@connect((state) => ({
  reembedParams: state.reembedParameters,
  annoMatrix: state.annoMatrix
}))
class DimredPanel extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      cfshown: false,
      gfshown: false,
      hvgshown: false,
      samshown: false,
      trshown: false,
    };
  }

  render() {
    const {
      cfshown, gfshown, hvgshown, samshown, trshown
    } = this.state;
    const { reembedParams, annoMatrix, dispatch, embName, onChange } = this.props;
    return (
      <div>
      <DefaultsButton dispatch={dispatch}/>  
      <div
      style={{
        paddingBottom: "10px",
        paddingTop: "10px"
      }}>
      <InputGroup
          id="emb-name-input"
          placeholder="New embedding name..."
          onChange={onChange}
          value={embName}
      />
      </div>
      <AnchorButton
        onClick={() => {
          this.setState({ 
            hvgshown: false,
            cfshown: false,
            gfshown: false,
            trshown: !this.state.trshown,
            samshown: false,
          });
        }}
        text={`Transformation`}
        fill outlined
        rightIcon={trshown ? "chevron-down" : "chevron-right"} small
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
            <ParameterInput
              label="Scale data?"
              param="scaleData"
            />
          </ControlGroup> 
        </Collapse>  
      </div>        
      <AnchorButton
        onClick={() => {
          this.setState({ 
            hvgshown: false,
            cfshown: !this.state.cfshown,
            gfshown: false,
            samshown: false,
            trshown: false,
          });
        }}
        text="PCA"
        fill outlined
        rightIcon={cfshown ? "chevron-down" : "chevron-right"} small
      />                    
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={cfshown}>
          <ControlGroup fill={true} vertical={false}>
            <ParameterInput
            min={1}
            label="n_comps"
            param="numPCs"
            />
            <div style={{"margin":"auto 0"}}>
              <ParameterInput
              label="svd_solver"
              param="pcaSolver"
              options={["arpack","randomized","auto","lopcg"]}
              />   
            </div>             
          </ControlGroup>
        </Collapse>
      </div>     
      <AnchorButton
        onClick={() => {
          this.setState({ 
            hvgshown: false,
            gfshown: !this.state.gfshown,
            cfshown: false,
            samshown: false,
            trshown: false,
          });
        }}
        text={`Neighbors`}
        fill outlined
        rightIcon={gfshown ? "chevron-down" : "chevron-right"} small
      />   
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={gfshown}>
          <ControlGroup fill={true} vertical={false}>
            <ParameterInput
              min={1}
              label="n_neighbors"
              param="neighborsKnn"
              />  
            <div style={{"margin":"auto 0"}}>
              <ParameterInput
                label="metric"
                param="distanceMetric"
                options={["euclidean","correlation","cosine"]}
                />
            </div>
            <div style={{"margin":"auto 0"}}>
              <ParameterInput
                label="method"
                param="neighborsMethod"
                options={["umap","gauss","rapids"]}
                />   
              </div>                                    
            </ControlGroup>         
        </Collapse>
      </div>
      <AnchorButton
        onClick={() => {
          this.setState({ 
            hvgshown: false,
            cfshown: false,
            gfshown: false,
            samshown: !this.state.samshown,
            trshown: false,
          });
        }}
        text={`SAM`}
        fill outlined
        rightIcon={samshown ? "chevron-down" : "chevron-right"} small
        disabled = {!reembedParams.doSAM}
      />   
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={samshown}>    
          <ControlGroup fill={true} vertical={false}>
            <ParameterInput
              label="num_norm_avg"
              param="nnaSAM"
              min={1}
              max={annoMatrix.schema.nVar}
            />     
            <div style={{"margin":"auto 0"}}>
              <ParameterInput
                label="num_norm_avg"
                param="weightModeSAM"
                options={["dispersion","rms","combined"]}
              />    
            </div> 
          </ControlGroup>                 
        </Collapse>  
      </div>         
      <AnchorButton
        onClick={() => {
          this.setState({ 
            hvgshown: !this.state.hvgshown,
            cfshown: false,
            gfshown: false,
            samshown: false,
          });
        }}
        text={`UMAP`}
        fill outlined
        rightIcon={hvgshown ? "chevron-down" : "chevron-right"} small
      />   
      <div style={{"paddingLeft":"10px"}}>
        <Collapse isOpen={hvgshown}>      
          <ParameterInput
            min={0.01}
            label="min_dist"
            param="umapMinDist"
          />            
        </Collapse>  
      </div>                  
    </div>
    );
  }
}

export default DimredPanel;
