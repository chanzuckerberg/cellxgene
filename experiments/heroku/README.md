# cellxgene cloud deployment with Heroku

## Quickstart

Clicking on the following button will forward you to Heroku to begin the deployment process:

<a href="https://heroku.com/deploy?template=https://github.com/chanzuckerberg/cellxgene/tree/heroku">
  <img src="https://www.herokucdn.com/deploy/button.svg" alt="Deploy">
</a>

If not already logged in to Heroku, there you will be prompted to log in or sign up for an account.

Once logged in you will be sent to the setup page. Here you can set some of the basic settings for the app:

#### Default settings

- `App name`: the unique name for your deployment
  - This will also serve as the default URL (e.g. https://cellxgene.herokapp.com/)
- `App owner`: Who will own this app. Either you personally or an organization/team
- `Region`: Location of the server where the app will be deployed (EU or US)

#### Configuration

- `DATASET`: A _publicly_ accessible URL pointing to a .h5ad file to view
  - This defaults to pbm3k.h5ad

After filling out the settings and pressing the `Deploy app` button Heroku will begin building your deployment. This process will take a few minutes, but once completed you will have a personal free hosted version of cell×gene!

## What is Heroku

Heroku is a service provided by Salesforce which offers a quick and easy substitute for hosting applications on the cloud.

A Heroku deployment of cell×gene means that the app is not running on your local machine. Instead, the app is installed, configured, and ran on the Heroku servers (read: cloud).

Heroku is one of many options available for hosting instances of cellxgene on the web.
Some other options include: Amazon Web Services, Google Cloud Platform, Digital Ocean, and Microsoft Azure

## Why should I use Heroku to deploy cellxgene?

Because cellxgene currently heavily relies on its Python backend for providing the client with the necessary data and tooling, it is currently not possible to host cellxgene as a static webpage via Netlify, AWS, GitHub Pages, etc.

What Heroku enables is a quick, non-technical method of setting up a cell×gene instance. No command line knowledge needed. This also allows machines anywhere to access the instance, so sharing a visualized dataset is as simple as sharing a link.

This is a good option if you want to quickly deploy an instance of cellxgene to the web. Heroku deployments are free for small datasets up to around 500MBs in size. See below regarding larger datasets.

## When should I not deploy with Heroku?

- The default free dyno offered by Heroku is limited in memory to 512 MBs
  - The amount of memory needed for the dyno is roughly the same size as the h5ad file
  - Heroku offers tiered paid dynos. More can be found [here](https://www.heroku.com/pricing)
  - Note that this can get _very_ expensive for larger datasets
- On the free dyno, after 30 minutes of inactivity, Heroku will put your app into a hibernation mode. Requiring it to boot up on the next access
- Having multiple simultaneous users requires more memory. This means that the free container size is easily overwhelmed by multiple users, even with small datasets; this can be addressed by purchasing a larger container size.
- For the dataset to be downloaded onto the distro, your dataset must be hosted on a publicly accessible URL
