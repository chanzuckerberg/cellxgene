## Creating PR
1. Name [username]/branchname 
    1. Branch name should be all lowercase
    2. Words separated by “-”
2. Code should address only one issue ideally, make a separate PR for each task
3. Description 
    1. Clear explanation of issues solved
    2. Describe why and how when appropriate
    3. Call out specific areas you want extra attention in review (optional)
    4. If your PR requires more than one reviewer tag those people in the description or comments and let them know that you specifically require them
    4. Ensure that the PR updates tests and documentation and adds tests where appropriate
4. Fixes issues
    1. Use github’s issue keywords when PR is addressing an issue https://help.github.com/articles/closing-issues-using-keywords/
5. Tags (add at beginning of title)
    1. [EASY] - small non-controversial change, easy to review
    2. [DO NOT MERGE] - PR is in progress, do not merge changes

## Review
1. Assign at least one reviewer to submitted PRs.   Reviewers should be selected based on expertise in areas affected by the PR (eg, web UI: Colin), and should include Comp Bio and PM as needed.
2. Reviewers should approve or request changes (not just comment) and put general and line level comments where appropriate
3. As a PR submitter respond to all comments (eg, comment, commit a change, etc)
4. External PRs
    1. For external PRs or PRs not from our core team, core team should assign a reviewer and make initial contact within 1 business day
    2. Build code on local environment and run smoke tests 

## Required to Merge
1. Travis CI Build passing
2. At least one reviewer approved
    1. Exceptions: 
        1. Release PRs where version is just bumped should not need review
        2. Complex PRs which touch multiple parts of the codebase should have reviews from all relevant parties
3. License and Security checks (SNYK) passing. If their server is down and you didn’t add any new external npm or python packages, merge is OK

## Merging
1. Use "squash and merge" option when merging
2. If you resolved conflicts, wait until the build passes to merge
