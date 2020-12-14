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
        if (retryCount === 0) {
          fetchAndHandle();
          // delay amount (3.5, 17.5, 87.5, 437.5 seconds)
        } else timeoutId = setTimeout(fetchAndHandle, 700 * 5 ** retryCount);
      }
      async function fetchAndHandle() {
        try {
          promise.response = await mockFetch(...fetchArgs);

          if (promise.response.ok || retryCount === retryTimes) {
            resolve(promise.response);
          }
          console.log("success", fetchArgs[2]);
        } catch (error) {
          console.log("error", fetchArgs[2]);
          if (retryCount === retryTimes) {
            reject(error);
          }
        }
        retryCount += 1;

        delayFetch();
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
    const response = { ok: false };
    if (rand > 5) {
      if (rand % 1 === 1) reject(new Error("fetch error"));
    } else response.ok = true;
    resolve(response);
  });
}
