export default function createBackoffFetch(retryTimes = 3) {
  let prevHash = null;
  let retryCount = 0;
  let promise = null;
  let timeoutId = null;

  return (hash, ...fetchArgs) => {
    if (prevHash !== hash) {
      prevHash = hash;
      retryCount = 0;

      if (timeoutId) {
        clearTimeout(timeoutId);
        timeoutId = null;
      }

      promise = backoffFetch(fetchArgs);
    }

    // This isn't an initial call, return the existing promise
    return promise;
  };

  function backoffFetch(fetchArgs) {
    return new Promise(promiseExecuter);

    function promiseExecuter(resolve, reject) {
      delayFetch();

      function delayFetch() {
        timeoutId = setTimeout(async () => {
          try {
            const response = await fetch(...fetchArgs);

            if (!response.ok) {
              const { status, statusText } = response;

              throw Error(`${status}: ${statusText}`);
            }

            resolve(response);
          } catch (error) {
            if (retryCount === retryTimes) {
              reject(error);
            }

            retryCount += 1;

            delayFetch();
          }
          // delay amount (3.5, 17.5, 87.5, 437.5 seconds)
        }, retryCount * 700 * 5 ** retryCount);
      }
    }
  }
}
