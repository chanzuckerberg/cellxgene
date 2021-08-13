export default function fromEntries<T = unknown>(
  arr: [string | number, T][]
): { [key: string]: T } {
  /*
	Similar to Object.fromEntries, but only handles array.
	This could be replaced with the standard function once it
	is widely available.   As of 3/20/2019, it has not yet
	been released in the Chrome stable channel.
	*/
  const obj: { [key: string]: T } = {};

  for (let i = 0; i < arr.length; i += 1) {
    obj[arr[i][0]] = arr[i][1];
  }

  return obj;
}
