# Heroku Deployment for cell×gene

## Quick Start

Clicking on the following button will forward you to heroku to begin the deployment process

<a href="https://heroku.com/deploy?template=https://github.com/chanzuckerberg/cellxgene/tree/heroku">
  <img src="https://www.herokucdn.com/deploy/button.svg" alt="Deploy">
</a>

If not already logged in, there you will be prompted to login or sign up for an account.

Once logged in you will be sent to the setup page. Here you can set some of the basic setting for the app:

#### Default Settngs

- `App name`: the unique name name for your deployment.
  - This will also serve as the default url (e.g. https://cellxgene.herokapp.com/)
- `App owner`: Who will own this app. Either you personally or an organiziation/team
- `Region`: Location of server where app will be deployed (EU or US)

#### Config Vars

- `DATASET`: A _publicly_ accessable URL pointing to a .h5ad file to view
  - This is defaulted to pbm3k.h5ad

After filling out the settings and pressing the `Deploy app` button heroku will begin building your deployment. This process will take a few minutes, but once completed you will have a personal free hosted version of cell×gene!

## Caveats

- The default free dyno offered by Heroku is limited in memory to 512 MBs
  - Heroku offers tiered paid dynos. More can be found [here](https://www.heroku.com/pricing
  - The amount of memory needed for the dyno is equal to the the size of the h5ad file
- Inorder for the dataset to be downloaded during app setup, your dataset must be hosted on a publicly accessable url
