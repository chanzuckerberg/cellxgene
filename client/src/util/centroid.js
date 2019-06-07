/*
  Centroid coordinate calculation


*/

export default function calcCentroid(
  world,
  annoName,
  annoValue,
  layoutDimNames
) {
  const annoArray = world.obsAnnotations.col(annoName).asArray();
  const layoutXArray = world.obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = world.obsLayout.col(layoutDimNames[1]).asArray();

  const centroidData = annoArray.reduce(
    (data, val, i) => {
      if (val === annoValue) {
        data.x += layoutXArray[i];
        data.y += layoutYArray[i];
        data.count += 1;
      }
      return data;
    },
    { x: 0, y: 0, count: 0 }
  );

  const centroid = [];

  centroid[0] = centroidData.x / centroidData.count;
  centroid[1] = centroidData.y / centroidData.count;
  centroid[2] = centroidData.count / world.nObs;
  

  return centroid;

  // Optimization from bruce, cut down on object creation/deletion
  // const {x,y,count } = centroidData;
  // x = x/count;
  // y = y/count;
  // return { x, y, count };

  // let x = 0;
  // let y = 0;
  // let count = 0;
  // for (let i =0, l = annoArray.length ; i < l; i += 1) {
  //   if (...) {
  //     x += layoutXarray[i];
  //   }

  // }
}
