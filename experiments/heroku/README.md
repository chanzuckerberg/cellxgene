# Heroku Deployment for cell×gene

## Quick Start

Clicking on the following button will forward you to heroku to begin the deployment process

<a href="https://heroku.com/deploy?template=https://github.com/chanzuckerberg/cellxgene/tree/heroku">
  <img src="https://www.herokucdn.com/deploy/button.svg" alt="Deploy">
</a>

If not already logged in to heroku, there you will be prompted to login or sign up for an account.

Once logged in you will be sent to the setup page. Here you can set some of the basic settings for the app:

#### Default Settngs

- `App name`: the unique name name for your deployment.
  - This will also serve as the default url (e.g. https://cellxgene.herokapp.com/)
- `App owner`: Who will own this app. Either you personally or an organiziation/team
- `Region`: Location of server where app will be deployed (EU or US)

#### Config Vars

- `DATASET`: A _publicly_ accessable URL pointing to a .h5ad file to view
  - This is defaulted to pbm3k.h5ad

After filling out the settings and pressing the `Deploy app` button heroku will begin building your deployment. This process will take a few minutes, but once completed you will have a personal free hosted version of cell×gene!

## What Is Heroku and Why Should I Use It?

> Heroku is a cloud platform that lets companies build, deliver, monitor and scale apps — we're the fastest way to go from idea to URL, bypassing all those infrastructure headaches.

A Heroku deployment of cell×gene means that the app is not running on your local machine. Instead, the app is installed, configured, and ran on the Heroku servers (read: cloud).

What this enables is a quick, non-technical method of setting up a cell×gene instance. No command line knowledge needed. This also allows machines anywhere to access the instance, so sharing a visualized dataset is as simple as sharing a link.

If you want to quickly deploy an instance of cellxgene and, if not using a small dataset, okay with paying for a tiered dyno, this is a viable solution.

## Why Shouldn't I Use It?

- The default free dyno offered by Heroku is limited in memory to 512 MBs
  - Heroku offers tiered paid dynos. More can be found [here](https://www.heroku.com/pricing)
  - The amount of memory needed for the dyno is equal to the the size of the h5ad file
  - This can get _very_ expensive for larger datasets
- Inorder for the dataset to be downloaded during app setup, your dataset must be hosted on a publicly accessable url
- The server is still easily overwhelmed by many concurrent users
  - This can be combatted by paying for a higher tiered dyno
