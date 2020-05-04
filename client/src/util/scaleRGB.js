// jshint esversion: 6
export default (input) => {
  const outputMax = 1;
  const outputMin = 0;

  const inputMax = 255;
  const inputMin = 0;

  const percent = (input - inputMin) / (inputMax - inputMin);
  return percent * (outputMax - outputMin) + outputMin;
};
