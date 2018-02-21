import _ from "lodash";

import {
  margin,
  width,
  height,
  createDimensions,
} from "./util";

const drawScatterplotCanvas = (
  context,
  xScale,
  yScale,
  currentCellSelection,
  opacityForDeselectedCells,
  expression,
  scatterplotXXaccessor,
  scatterplotYYaccessor
) => {

  /* clear canvas */
  context.clearRect(0, 0, width, height);
  
  const _currentCellSelectionMap = _.keyBy(currentCellSelection, "CellName"); /* move me to the reducer */

  expression.data.cells.forEach((cell, i) => {
    
    /* 
      this is necessary until we are no longer getting expression for all cells, but only for 'world'
      ...which will mean refetching when we regraph, or 'go back up to all cells'
    */
      
    if (!_currentCellSelectionMap[cell.cellname]) { return } 
    
    context.beginPath();
    /* context.arc(x,y,r,sAngle,eAngle,counterclockwise); */
    context.arc(
      xScale(cell.e[expression.data.genes.indexOf(scatterplotXXaccessor)]),            /* x */
      yScale(cell.e[expression.data.genes.indexOf(scatterplotYYaccessor)]),            /* y */
      _currentCellSelectionMap[cell.cellname]["__selected__"] ? 3 : 1.5,                   /* r */
      0,                  /* sAngle */
      2 * Math.PI         /* eAngle */
    );
    
    context.fillStyle = _currentCellSelectionMap[cell.cellname]["__color__"]

    if (_currentCellSelectionMap[cell.cellname]["__selected__"]) {
      context.globalAlpha = 1;
    } else {
      context.globalAlpha = opacityForDeselectedCells;
    }

    context.fill();
  });
}

export default _.debounce(drawScatterplotCanvas, 50)