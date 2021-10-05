import React from "react";
import { connect } from "react-redux";
import { useAsync } from "react-async";
import {
  Button,
  ButtonGroup,
  H4,
  Popover,
  Position,
  Radio,
  RadioGroup,
  AnchorButton,
  Tooltip,
} from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import { getDiscreteCellEmbeddingRowIndex } from "../../util/stateManager/viewStackHelpers";
import _ from "lodash";

@connect((state) => {
  return {
    layoutChoice: state.layoutChoice, // TODO: really should clean up naming, s/layout/embedding/g
    schema: state.annoMatrix?.schema,
    annoMatrix: state.annoMatrix,
    crossfilter: state.obsCrossfilter,
  };
})
class Embedding extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {};
  }

  handleLayoutChoiceChange = (e) => {
    const { dispatch, layoutChoice } = this.props;
    if (layoutChoice.available.includes(e.currentTarget.value)) {
      dispatch(actions.layoutChoiceAction(e.currentTarget.value));
    }
  };
  handleDeleteEmbedding = (e,val) => {
    const { dispatch, annoMatrix, layoutChoice } = this.props;
    const { available } = layoutChoice;
    const toDelete = [val]
    available.forEach((item) => {
      if (item.includes(`${val};;`)){
        toDelete.push(item)
      }
    });
    let newAnnoMatrix;
    toDelete.forEach((item) => {
      dispatch({type: "reembed: delete reembedding", embName: item})
      newAnnoMatrix = annoMatrix.dropObsmLayout(val);
      dispatch({type: "", annoMatrix: newAnnoMatrix})
    })
    dispatch(actions.requestDeleteEmbedding(toDelete))
  }
  render() {
    const { layoutChoice, schema, crossfilter } = this.props;
    const { annoMatrix } = crossfilter;
    return (
      <ButtonGroup
        style={{
          position: "absolute",
          display: "inherit",
          left: 8,
          bottom: 8,
          zIndex: 9999,
        }}
      >
        <Popover
          target={
            <Tooltip
              content="Select embedding for visualization"
              position="top"
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
            >
              <Button
                type="button"
                data-testid="layout-choice"
                icon="heatmap"
                // minimal
                id="embedding"
                style={{
                  cursor: "pointer",
                }}
              >
                {layoutChoice?.current.split(";;").at(-1)}: {crossfilter.countSelected()} out of{" "}
                {crossfilter.size()} cells
              </Button>
            </Tooltip>
          }
          // minimal /* removes arrow */
          position={Position.TOP_LEFT}
          content={
            <div
              style={{
                display: "flex",
                justifyContent: "flex-start",
                alignItems: "flex-start",
                flexDirection: "column",
                padding: 10,
                width: 400,
              }}
            >
              <H4>Embedding Choice</H4>
              <p style={{ fontStyle: "italic" }}>
                There are {schema?.dataframe?.nObs} cells in the entire dataset.
              </p>
              <EmbeddingChoices
                onChange={this.handleLayoutChoiceChange}
                annoMatrix={annoMatrix}
                layoutChoice={layoutChoice}
                onDeleteEmbedding={this.handleDeleteEmbedding}
              />
            </div>
          }
        />
      </ButtonGroup>
    );
  }
}

export default Embedding;

const loadAllEmbeddingCounts = async ({ annoMatrix, available }) => {
  const embeddings = await Promise.all(
    available.map((name) => annoMatrix.base().fetch("emb", name))
  );
  return available.map((name, idx) => ({
    embeddingName: name,
    embedding: embeddings[idx],
    discreteCellIndex: getDiscreteCellEmbeddingRowIndex(embeddings[idx]),
  }));
};

const EmbeddingChoices = ({ onChange, annoMatrix, layoutChoice, onDeleteEmbedding }) => {
  const { available } = layoutChoice;
  const { data, error, isPending } = useAsync({
    promiseFn: loadAllEmbeddingCounts,
    annoMatrix,
    available,
  });
  
  if (error) {
    /* log, as this is unexpected */
    console.error(error);
  }
  if (error || isPending) {
    /* still loading, or errored out - just omit counts (TODO: spinner?) */
    return (
      <div>
        {null}
      </div>
    );
  }
  if (data) {
    const name = layoutChoice.current;
    let parentName;
    const embName = name;
    if(name.includes(";;")){
      parentName = name.replace(`;;${name.split(";;").at(-1)}`,"")
    } else {
      parentName = "";
    }
    let sizeHintParent;
    let sizeHintCurrent;
    const x = data.map((summary) => {
      const { discreteCellIndex, embeddingName } = summary;
      const sizeHint = `${discreteCellIndex.size()} cells`;

      let queryParent;
      const queryName = embeddingName;
      if(embeddingName.includes(";;")){        
        queryParent = embeddingName.replace(`;;${embeddingName.split(";;").at(-1)}`,"")
      } else {
        queryParent = "";
      }      
      if (queryName === embName) {
        sizeHintCurrent = sizeHint;
      }       
      if (queryName === parentName) {
        sizeHintParent = sizeHint;
      }
      if((parentName === "" && queryParent === "" || queryParent === embName) && available.includes(queryName)){
        return (
         
            <Radio
              label={`${queryName.split(';;').at(-1)}: ${sizeHint}`}
              value={embeddingName}
              key={embeddingName}
              style={{
                display: "flex",
                verticalAlign: "middle",
              }}
              children={
                (queryName !== embName ? <AnchorButton
                  icon="small-cross"
                  minimal
                  style={{
                    cursor: "pointer",
                    marginLeft: "auto",
                    marginTop: "-5px"
                  }}
                  onClick={(e) => onDeleteEmbedding(e,queryName)}
                /> : null)
              }
            />  
            
        // add `X` button which deletes embedding for all embeddings that are not present in original schema.
        // remove all embeddings that are its children 
        // to figure out: how to delete embeddings from schema.
        );
      } else {
        return null;
      }
    });
    if(parentName !== ""){
      x.unshift(
        <Radio
          label={`${embName.split(';;').at(-1)}: ${sizeHintCurrent}`}
          value={name}
          key={name}
        />        
      )
      x.unshift(
        <Radio
          label={`(Parent) ${parentName.split(';;').at(-1)}: ${sizeHintParent}`}
          value={parentName}
          key={parentName}
        />        
      )      
    }

    return (
      <RadioGroup onChange={onChange} selectedValue={layoutChoice.current}>
        {x}
      </RadioGroup>
    );
  }
  return null;
};
