import React from "react";
import { connect } from "react-redux";
import {
  Button,
  ButtonGroup,
  H4,
  Popover,
  Position,
  Radio,
  RadioGroup,
  AnchorButton,
  MenuItem,
  Menu,
  Tooltip,
} from "@blueprintjs/core";
import * as globals from "../../globals";
import actions from "../../actions";
import { getDiscreteCellEmbeddingRowIndex } from "../../util/stateManager/viewStackHelpers";
import _ from "lodash";
import AnnoDialog from "../annoDialog";
import LabelInput from "../labelInput";

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
    this.state = {newLayoutText: "", isEmbeddingExpanded: {"": true}};
  }
  
  handleChangeOrSelect = (name) => {
    this.setState({
      newLayoutText: name,
    });
  };
  closeAllChildren = (node,tree,expanded) => {
    const children = tree[node]?.children;    
    if (children){
      expanded[node] = false;
      for (const child of children){
        this.closeAllChildren(child,tree,expanded)
      }
    }
  }
  handleEmbeddingExpansionChange = (e,node,val,tree) => {
    const { isEmbeddingExpanded: newExpanded } = this.state;
    if (val) {
      this.closeAllChildren(node,tree,newExpanded)
      this.setState({
        isEmbeddingExpanded: {...newExpanded}
      })      
    } else {
      this.setState({
        isEmbeddingExpanded: {...newExpanded, [node]: true}
      })
    }
    if(e){
      e.preventDefault()
    }
  }
  activateEditLayoutMode = (e, embeddingName) => {
    const { dispatch } = this.props;
    this.setState({
      newLayoutText: embeddingName.split(';;').at(-1)
    })
    dispatch({
      type: "reembed: activate layout edit mode",
      data: embeddingName,
    });
  };
  disableEditLayoutMode = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "reembed: deactivate layout edit mode",
    });
  };
  handleEditLayout = (e) => {
    const { dispatch, layoutChoice } = this.props;
    const { newLayoutText, isEmbeddingExpanded } = this.state
    const { available } = layoutChoice;
    const toRename = [layoutChoice.layoutNameBeingEdited]
    available.forEach((item) => {
      if (item.includes(`${layoutChoice.layoutNameBeingEdited};;`)){
        toRename.push(item)
      }
    });
    const oldName = layoutChoice.layoutNameBeingEdited.split(';;').at(-1);
    const newName = newLayoutText;
    
    toRename.forEach((item) => {
      dispatch({type: "reembed: rename reembedding", embName: item, newName: item.replace(oldName,newName)})
      const val = isEmbeddingExpanded?.[item] ?? false;
      this.setState({
        isEmbeddingExpanded: {...isEmbeddingExpanded, [item]: _, [item.replace(oldName,newName)]: val}
      })
    })
    if (oldName !== newName) {
      dispatch(actions.requestRenameEmbedding(toRename,oldName,newName))
    }
    dispatch({
      type: "reembed: deactivate layout edit mode",
    });       
  }
  handleLayoutChoiceChange = (e) => {
    const { dispatch, layoutChoice } = this.props;
    const { isEmbeddingExpanded } = this.state
    if (layoutChoice.available.includes(e.currentTarget.value) && !layoutChoice.isEditingLayoutName) {
      this.setState({
        isEmbeddingExpanded: {...isEmbeddingExpanded, [e.currentTarget.value]: true}
      })
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
    const { newLayoutText, isEmbeddingExpanded } = this.state;
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
                activateEditLayoutMode={this.activateEditLayoutMode}
                isEmbeddingExpanded={isEmbeddingExpanded}
                handleEmbeddingExpansionChange={this.handleEmbeddingExpansionChange}
              />
              <AnnoDialog
                isActive={
                  layoutChoice.isEditingLayoutName
                }
                inputProps={{
                  "data-testid": `edit-layout-name-dialog`,
                }}
                primaryButtonProps={{
                  "data-testid": `submit-layout-edit`,
                }}
                title="Edit layout name"
                instruction={"Choose a new layout name"}
                cancelTooltipContent="Close this dialog without editing this layout."
                primaryButtonText="Edit layout name"
                text={layoutChoice.layoutNameBeingEdited}
                handleSubmit={this.handleEditLayout}
                handleCancel={this.disableEditLayoutMode}
                annoInput={
                  <LabelInput
                    label={newLayoutText}
                    inputProps={{
                      "data-testid": `edit-layout-name-text`,
                      leftIcon: "tag",
                      intent: "none",
                      autoFocus: true,
                    }}
                    onChange={this.handleChangeOrSelect}
                    onSelect={this.handleChangeOrSelect}                    
                    newLabelMessage="New layout name"
                  />
                }
              />              
            </div>
          }
        />
      </ButtonGroup>
    );
  }
}

export default Embedding;

const loadAllEmbeddingCounts = async (annoMatrix, available) => {
  const embeddings = await Promise.all(
    available.map((name) => annoMatrix.base().fetch("emb", name))
  );
  try {
    return available.map((name, idx) => ({
      embeddingName: name,
      embedding: embeddings[idx],
      discreteCellIndex: getDiscreteCellEmbeddingRowIndex(embeddings[idx]),
    }));
  } catch {
    // nothing happens
  }
  
};

/*
below function will generate an array with full tree of embeddings.
but i only want to show PARENT+children, selecting parent should switch to parent's parent. if no parent, padding = 0

so i can use current indented embedding tree, but the node I give it is going to change based on user's selection.

*/

const IndentedEmbeddingTree = (node,roots,tree,padding, els, currView, onDeleteEmbedding, activateEditLayoutMode, isEmbeddingExpanded, handleEmbeddingExpansionChange) => {
  const children = tree[node]?.children;
  els.push(
    (isEmbeddingExpanded?.[tree[node].parent] ?? tree[tree[node].parent].expandedByDefault) ? <Radio
      label={`${node.split(';;').at(-1)}: ${tree[node].sizeHint}`}
      value={node}
      key={node}
      style={{
        display: "flex",
        verticalAlign: "middle",
        paddingLeft: `${padding+26}px`
      }}
      children={
      <div style={{
        paddingLeft: "5px",
      }}>
      {children && !(tree[node]?.disable ?? false) ? 
      <AnchorButton
        icon={isEmbeddingExpanded?.[node] ?? tree[node].expandedByDefault ? "chevron-down" : "chevron-right"}
        data-testid={`${node}:expand-embeddings`}
        onClick={(e) => handleEmbeddingExpansionChange(e,node,isEmbeddingExpanded?.[node] ?? tree[node].expandedByDefault, tree)}
        minimal
        style={{
          cursor: "pointer",
          marginLeft: "auto",
          marginTop: "-5px"
        }}                    
        /> : null}

        {node !== "root" ?
      <AnchorButton
          icon="more"
          data-testid={`${node}:edit-layout-mode`}
          onClick={(e) => activateEditLayoutMode(e,node)}
          minimal
          style={{
            cursor: "pointer",
            marginLeft: "auto",
            marginTop: "-5px"
          }}                    
        /> : null}
        {node !== currView && node !== "root"  && !roots.includes(node) ?

        <AnchorButton
          icon="small-cross"
          minimal
          style={{
            cursor: "pointer",
            marginLeft: "auto",
            marginTop: "-5px"
          }}
          onClick={(e) => onDeleteEmbedding(e,node)}
        /> : null}  
      </div> 
      }
    />  : null
  )  
  if (children){
    for (const child of children){
      IndentedEmbeddingTree(child,roots,tree,padding+26, els, currView, onDeleteEmbedding, activateEditLayoutMode, isEmbeddingExpanded, handleEmbeddingExpansionChange)
    }
  }
}

const EmbeddingChoices = ({ onChange, annoMatrix, layoutChoice, onDeleteEmbedding, activateEditLayoutMode, isEmbeddingExpanded, handleEmbeddingExpansionChange }) => {
  const [ data, setData ] = React.useState(null)
  const [ renderedEmbeddingTree, setRenderedEmbeddingTree ] = React.useState(null)
  React.useEffect(() => {
    const { available } = layoutChoice;
    loadAllEmbeddingCounts(annoMatrix,available).then((res)=>{
      setData(res)
      if (res) {
        const name = layoutChoice.current;
        let parentName;
        if(name.includes(";;")){
          parentName = name.replace(`;;${name.split(";;").at(-1)}`,"")
        } else {
          parentName = "";
        }
    
        const embeddingTree = {}
        res.map((summary) => {
          const { discreteCellIndex, embeddingName: queryName } = summary;
          let queryParent;
          if(queryName.includes(";;")){        
            queryParent = queryName.replace(`;;${queryName.split(";;").at(-1)}`,"")
          } else {
            queryParent = "";
          }
          
          const sizeHint = `${discreteCellIndex.size()} cells`;
          // add queryName to children of queryParent
          if (embeddingTree?.[queryParent]?.children) { //if children exists on queryParent
            embeddingTree[queryParent].children.push(queryName)
          } else if (embeddingTree?.[queryParent]) { // if anything exists on queryParent
            embeddingTree[queryParent] = {...embeddingTree[queryParent], children: [queryName]}
          } else { // create new entry for queryParent
            embeddingTree[queryParent] = {children: [queryName]}
          }

          const expandedByDefault = (queryParent==="" || queryName === name);

          if (embeddingTree?.[queryName]){ // queryName exists in embeddingTree
            embeddingTree[queryName] = {...embeddingTree[queryName], sizeHint: sizeHint, expandedByDefault: expandedByDefault, parent: queryParent}
          } else {
            embeddingTree[queryName] = {sizeHint: sizeHint, expandedByDefault: expandedByDefault, parent: queryParent}
          }
        });      
        const els = []
        let iterable;
        let roots;
        if (parentName === ""){
          iterable = embeddingTree[""].children      
        } else {
          let currNode = parentName;
          let iterate = true;
          roots = [currNode]
          embeddingTree[currNode].expandedByDefault = true;
          embeddingTree[currNode].disable = true;
          while (iterate){
            if (embeddingTree[currNode].parent === ""){
              iterate=false;
            } else {
              currNode = embeddingTree[currNode].parent
              embeddingTree[currNode].expandedByDefault = true;
              embeddingTree[currNode].disable = true;
              roots.push(currNode)
            }
          }
          iterable = [currNode]
        }
        
        for (const c of iterable){
          IndentedEmbeddingTree(c,roots??iterable.filter((item)=>item!==c),embeddingTree,0, els, name, onDeleteEmbedding, activateEditLayoutMode, isEmbeddingExpanded, handleEmbeddingExpansionChange)    
        }
        setRenderedEmbeddingTree(
          <RadioGroup onChange={onChange} selectedValue={layoutChoice.current}>
            {els}
          </RadioGroup>
        );
      }
    })
  }, [annoMatrix, layoutChoice, isEmbeddingExpanded]);

  if (!data) {
    /* still loading, or errored out - just omit counts (TODO: spinner?) */
    return (
      <div>
        {null}
      </div>
    );
  }
  return renderedEmbeddingTree;

}






