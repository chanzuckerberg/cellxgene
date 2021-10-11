import React from "react";
import { connect } from "react-redux";
import { Button, ControlGroup } from "@blueprintjs/core";
import ParameterInput from "../menubar/parameterinput";
import GeneSet from "./geneSet";
import { GenesetHotkeys } from "../hotkeys";
import actions from "../../actions";
import CreateGenesetDialogue from "./menus/createGenesetDialogue";
import * as globals from "../../globals";
import { AnnoMatrixLoader } from "../../annoMatrix";
@connect((state) => {
  return {
    genesets: state.genesets.genesets,
    colorAccessor: state.colors.colorAccessor,
    annoMatrix: state.annoMatrix,
    reembedParams: state.reembedParameters
  };
})
class GeneExpression extends React.Component {
  renderGeneSets = () => {
    const sets = [];
    const { genesets } = this.props;
    for (const [name, geneset] of genesets) {
      sets.push(
        <GeneSet
          key={name}
          setGenes={Array.from(geneset.genes.keys())}
          setGenesWithDescriptions={geneset.genes}
          setName={name}
          genesetDescription={geneset.genesetDescription}
        />
      );
    }
    return sets;
  };

  componentDidUpdate(prevProps) {
    const { dispatch, reembedParams, annoMatrix } = this.props;
    if(prevProps.reembedParams.dataLayerExpr !== reembedParams.dataLayerExpr){
      // Trigger new data layer.
      dispatch(actions.requestDataLayerChange(reembedParams.dataLayerExpr)).then(()=>{
        const baseDataUrl = `${globals.API.prefix}${globals.API.version}`;
        const annoMatrixNew = new AnnoMatrixLoader(baseDataUrl, annoMatrix.schema);
        dispatch({
          type: "",
          annoMatrix: annoMatrixNew
        });      
      })
    }
  }

  handleActivateCreateGenesetMode = () => {
    const { dispatch } = this.props;
    dispatch({ type: "geneset: activate add new geneset mode" });
  };

  render() {
    const { dispatch, genesets, colorAccessor, annoMatrix, reembedParams } = this.props;
    return (
      <div>
        <GenesetHotkeys
          dispatch={dispatch}
          genesets={genesets}
          colorAccessor={colorAccessor}
        />
        <div>
          <div style={{ display: "flex", marginBottom: 10, justifyContent: "space-between", position: "relative", top: -2 }}>
            <Button
              data-testid="open-create-geneset-dialog"
              onClick={this.handleActivateCreateGenesetMode}
              intent="primary"
            >
              Create new <strong>gene set</strong>
            </Button>
            <ParameterInput
              label="Data layer"
              param="dataLayerExpr"
              options={annoMatrix.schema.layers}
            />                   
                              
          </div>
          <CreateGenesetDialogue />
        </div>
        <div>{this.renderGeneSets()}</div>
      </div>
    );
  }
}

export default GeneExpression;
