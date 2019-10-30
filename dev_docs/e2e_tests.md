## End to end testing thesis
End to end (e2e) tests, when done correctly, can surface errors and regressions in the application before it is released to production. The tests are run against a version of the web app that is running the full stack. By checking the fully functional app, you are able to surface issues with logic and data at every level of the stack (you still have to figure out which level that error is occurring at). You are also able to test how features that are added affect the whole application.

The major problem most developers have with e2e tests are that they the tests can be flaky and thus not trusted. Failing tests wind up being disabled instead of fixed. To guard against that we need to establish principals about what features to test and have guidelines for writing tests to ensure that they are as robust as possible. By ensuring each test has a high value and is unlikely to fail when it shouldn't it increases the motivation to fix the test as opposed to deleting it.

### Choosing what to test:

* Common actions - Features and actions that are highly used will cause more impact if they break 
* High value - Features that may be rarely used, but cause a major business impact if they break 
* Depend on the whole stack and not just the FE - We want to test the backend and data as well as the FE. Features that only depend on the javascript code should be tested in unit tests instead.
* Feature is in a stable state - If the feature is in active development, it's likely that the test will have to be constantly updated. Wait until it is in a steady state before adding an e2e test.
* Feature is a likely target of regressions - If a feature has caused regressions in the past, it's good to add some monitoring to it.

### Guidelines for writing tests:

* Use data attributes as selectors instead of positional, classname, or id selectors. Add `data-testid` or `data-testclass` attributes to components that are inputs and outputs. The advantages are 1) You can move the component around and even change it's html element type and it will not break the test and 2) It is self documenting what the attribute is used by as opposed to a class name or id that could be used by CSS or something else.
* Use the helper functions in `puppeteerUtils.js` and ` for performing common actions. Puppeteer operates at a low level and these abstract it to a level of user actions.

For actions that don't meet these guidelines, you can create a test in feature.test.js that does not run on PRs so if the feature breaks it doesn't block anyone.

## Where/How
Where: client/\_\_tests\_\_/e2e/e2e.test.js

How: `npm run e2e`

Tests are written using the Jest testing framework. Browser automation is handled by puppeteer. Puppeteer provides an API to control Chrome over the DevTools Protocol. This allows us to write code that mimics user actions and check page state.

Tests are running on CI with every PR, they can also be run locally (see [developer guidelines](developer_guidelines.md)).  

The cellxgene instance is running with the pbmc3k dataset. Tests for additional datasets can be added, the data.js file contains the values for inputs and outputs expected for different datasets. Choose values that make sense for that dataset. 

## Setting up e2e tests locally
See [developer guidelines](developer_guidelines.md)

## How to write a new test
1. Ensure that your new test meets the guidelines above for high-value robust user actions to test.
2. Clearly define the series of user actions and outputs you want to test
3. Add testid (and/or testclass if element is not unique) data attributes to inputs and outputs you want to interact with.
4. Create a new test block possible nested in a describe block if you want to add more than one related test client/\_\_tests\_\_/e2e/e2e.test.js
5. Write the test code.
6. Run tests locally to ensure they pass.

**Debugging Tips** 
* Mystery bugs are often caused by a missing `async` or `await`. Double check every function and every line (twice!).
* Switch from to debug mode by running `export JEST_ENV='debug'` on the CLI before running the tests
* Put a javascript breakpoint in by adding this line `await jestPuppeteer.debug()` to the test at the point which you want to pause execution.
* Limit the tests that are run by adding `.only` after `describe` or `test` to limit execution to just that block of tests. (Remember to take the `.only` out before you commit!) 