export default function createBackoffFetch(retryTimes = 3) {
  let prevHash = null;
  let timeoutRef = null;
  let retryCount = 0;
  let response = null;
  let error = null;
  let promise = null;

  // This will only be called by itself or when we've seen a new matrix
  const backoffFetch = (fetchArgs) => {
    const promiseExecuter = (resolve, reject) => {
      retryCount += 1;
      // We don't want to timeout our first request
      if (retryCount === 1) {
        fetch(...fetchArgs)
          .then((resp) => {
            response = resp;
            if (resp.ok) {
              resolve(response);
            } else {
              promiseExecuter(resolve, reject);
            }
          })
          .catch((e) => {
            response = null;
            error = e;
            promiseExecuter(resolve, reject);
          });
      } else {
        // Timeout all other requests
        timeoutRef = setTimeout(async () => {
          // if (done) return;
          try {
            const resp = await fetch(...fetchArgs);
            response = resp;
            // Fetch status 2xx
            if (resp.ok) {
              resolve(response);
            } else if (retryCount - 1 === retryTimes) {
              if (error) reject(error);
              reject(response);
            } else {
              // Some non 2xx response code
              promiseExecuter(resolve, reject);
            }
          } catch (e) {
            // Fetch API error
            error = e;
            if (retryCount - 1 === retryTimes) reject(error);
            promiseExecuter(resolve, reject);
          }
          // delay amount (3.5, 17.5, 87.5, 437.5 seconds)
        }, 700 * 5 ** retryCount - 1);
      }
    };
    return new Promise(promiseExecuter);
  };

  return (hash, ...fetchArgs) => {
    if (prevHash !== hash) {
      clearTimeout(timeoutRef);
      prevHash = hash;
      timeoutRef = null;
      retryCount = 0;
      response = null;
      error = null;
      promise = backoffFetch(fetchArgs);
    }
    // This isn't an initial call, return the existing promise
    return promise;
  };
}
