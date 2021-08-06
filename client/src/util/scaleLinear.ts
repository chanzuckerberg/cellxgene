// Substitute for a d3 linear scale - less flexible, more performant.
// Returns a function which will scale a value.
//
// Example will scale [0,1] to [-1,1]
//     var myScale = scaleLinear([0, 1], [-1, 1]);
//     myScale(0) === -1
// this is is equivalent to d3.scaleLinear().domain([0,1]).range([-1,1])

export default (domain: Array<number>, range: Array<number>) => {
  const domainStart = domain[0];
  const scale = (range[1] - range[0]) / (domain[1] - domain[0]);
  const invScale = 1 / scale;
  const rangeStart = range[0];
  const f = (value: number) => (value - domainStart) * scale + rangeStart;

  // inverter
  f.invert = (value: number) => (value - rangeStart) * invScale + domainStart;

  return f;
};
