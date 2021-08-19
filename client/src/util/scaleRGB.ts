// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export default (input: any) => {
  const outputMax = 1;
  const outputMin = 0;

  const inputMax = 255;
  const inputMin = 0;

  const percent = (input - inputMin) / (inputMax - inputMin);
  return percent * (outputMax - outputMin) + outputMin;
};
