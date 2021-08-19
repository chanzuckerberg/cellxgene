export default function fromEntries(arr: any) {
  /*
	Similar to Object.fromEntries, but only handles array.
	This could be replaced with the standard fucnction once it
	is widely available.   As of 3/20/2019, it has not yet
	been released in the Chrome stable channel.
	*/
  const obj = {};
  for (let i = 0, l = arr.length; i < l; i += 1) {
    // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
    obj[arr[i][0]] = arr[i][1];
  }
  return obj;
}
