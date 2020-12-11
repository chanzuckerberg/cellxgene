const initReturnObject = {
  response: null,
  error: null,
  inProgress: false,
};
export default function createBackoffFetch(retryTimes = 3) {
  let prevHash = null;
  let timeoutRef = null;
  let retryCount = 0;
  let returnObject = initReturnObject;

  // This will only be called by itself or when we've seen a new matrix
  const backoffFetch = (fetchArgs) => {
    console.log(returnObject);
    // If we've expended all of our attempts
    if (retryCount >= retryTimes) {
      returnObject.inProgress = false;
      return returnObject;
    }
    // If we've marked the fetch as complete but haven't used up all of our attempts
    if (retryCount > 0 && !returnObject.inProgress) return returnObject;
    retryCount += 1;
    // DEBUG
    // DEBUG
    // DEBUG
    console.log("DELAY CALL anno", ...fetchArgs);
    console.log("DELAY CALL retryCount", retryCount);
    // We don't want to timeout our first request
    if (retryCount === 1) {
      returnObject.inProgress = true;
      fetch(...fetchArgs)
        .then((response) => {
          returnObject.response = response;
          if (response.ok) {
            returnObject.inProgress = false;
          }
        })
        .catch((e) => {
          returnObject.error = e;
        });
      return backoffFetch(fetchArgs);
    }
    // Timeout all other requests
    timeoutRef = setTimeout(async () => {
      try {
        const response = await fetch(...fetchArgs);
        returnObject.response = response;
        // Fetch status 2xx
        if (response.ok) {
          returnObject.inProgress = false;
        } else {
          // Some non 2xx response code
          return backoffFetch(fetchArgs);
        }
      } catch (e) {
        // Fetch API error
        returnObject.error = e;
        return backoffFetch(fetchArgs);
      }
      // Only way we get here is if there was no error - return return obj
      return returnObject;
      // delay amount (3.5, 17.5, 87.5, 437.5 seconds)
    }, 700 * 5 ** retryCount - 1);

    return returnObject;
  };

  return (hash, ...fetchArgs) => {
    // DEBUG
    // DEBUG
    // DEBUG
    console.log("-------prevAnno", prevHash);
    console.log("-------anno", hash);
    if (prevHash !== hash) {
      clearTimeout(timeoutRef);
      prevHash = hash;
      timeoutRef = null;
      retryCount = 0;
      returnObject = initReturnObject;
      backoffFetch(fetchArgs);
    }
    // This isn't an initial call, so just give a status check in the form of our return object
    return returnObject;
  };
}
