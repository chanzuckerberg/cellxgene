/*
  Centroid coordinate calculation


*/

export default function calcCentroid(world, annoName, annoValue, layoutDimNames) {
  const annoArray = world.obsAnnotations.col(annoName).asArray();
  const layoutXArray = world.obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = world.obsLayout.col(layoutDimNames[1]).asArray();

  const centroidData = annoArray.reduce((data, val, i) => {
    if(val === annoValue) {
      data.x += layoutXArray[i];
      data.y += layoutYArray[i];
      data.count += 1;
    }
    return data;
  }, { x: 0, y: 0, count: 0 });

  const centroid = { x: 0, y: 0 };

  centroid.x = centroidData.x / centroidData.count;
  centroid.y = centroidData.y / centroidData.count;

  return centroid;
}