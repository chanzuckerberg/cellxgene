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
          console.log("args:", fetchArgs[2], "attempt:", retryCount);
          try {
            const response = await mockFetch(...fetchArgs);

            if (!response.ok) {
              const { status, statusText } = response;

              throw Error(`${status}: ${statusText}`);
            }
            console.log("success", fetchArgs[2]);

            resolve(response);
          } catch (error) {
            console.log("error", fetchArgs[2]);
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

// DEBUG
// DEBUG
// DEBUG
function mockFetch() {
  return new Promise((resolve, reject) => {
    const rand = Math.floor(Math.random() * 10) + 1;
    const response = {};
    setTimeout(() => {
      if (rand > 5) {
        if (rand % 2) reject(new Error("fetch error"));
        response.ok = false;
      } else response.ok = true;
      resolve(response);
    }, rand * 500);
  });
}
