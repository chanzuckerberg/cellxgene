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

  const backoffFetch = (fetchArgs) => {
    if (retryCount >= retryTimes) {
      returnObject.inProgress = false;
      return returnObject;
    }
    if (retryCount > 0 && !returnObject.inProgress) return returnObject;
    retryCount += 1;
    // DEBUG
    // DEBUG
    // DEBUG
    console.log("DELAY CALL anno", ...fetchArgs);
    console.log("DELAY CALL retryCount", retryCount);
    if (retryCount === 1) {
      returnObject.inProgress = true;
      fetch(...fetchArgs)
        .then((response) => {
          returnObject.response = response;
          if (response.ok) {
            returnObject.inProgress = false;
            return returnObject;
          }
        })
        .catch((e) => {
          returnObject.error = e;
        });
      backoffFetch(fetchArgs);
    }
    timeoutRef = setTimeout(async () => {
      try {
        const response = await fetch(...fetchArgs);
        returnObject.response = response;
        if (response.ok) {
          returnObject.inProgress = false;
        } else {
          backoffFetch(fetchArgs);
        }
      } catch (e) {
        returnObject.error = e;
        backoffFetch(fetchArgs);
      }
    }, 700 * 5 ** retryCount - 1);

    return returnObject;
  };

  return (hash, ...fetchArgs) => {
    // DEBUG
    // DEBUG
    // DEBUG
    console.log("-------prevAnno", prevHash);
    console.log("-------anno", hash);
    if (!prevHash) {
      prevHash = hash;
    }
    if (prevHash !== hash) {
      clearTimeout(timeoutRef);
      prevHash = hash;
      timeoutRef = null;
      retryCount = 0;
      returnObject = initReturnObject;
      backoffFetch(fetchArgs);
    }
    return returnObject;
  };
}
