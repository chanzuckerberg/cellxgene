# Hosting cellxgene on the web

Cellxgene is intended to be used by researchers on their local machines. However, we recognize that sharing and exploring data on the web is important. We're exploring how we could better support this in the future, and [would welcome your input](https://github.com/chanzuckerberg/cellxgene/issues/875)!

In the meantime, you can see examples of how other groups have approached this in our [gallery](gallery). While we don't officially support web deployment, we've offered some guidance below on one way to deploy cellxgene to the web.

## General notes and cautions

Please consider the following when deploying cellxgene in any "hosted" environment, especially where access from the broader Internet is possible:

- Information security requires careful configuration of the host environment, including firewall, logging, etc. Please follow best practices.
- cellxgene includes features which may be inappropriate for a hosted deployment. You may wish to use the following command line option: `--disable-diffexp`.
- `cellxgene launch` currently uses Flask's development server, which is not recommended for hosted deployment (see the [Flask documentation](https://flask.palletsprojects.com/en/1.1.x/tutorial/deploy/#run-with-a-production-server))
- We have no testing or official support for deployments where multiple users are accessing the same cellxgene instance.
- Your cellxgene instance is likely to hang or crash if too many people access it at the same time, especially if they using functions that call the Python backend (such as differential expression, noted above).
- cellxgene only supports one instance per dataset

If you believe you have found a security-related issue with cellxgene, please report the issue immediately to <security@chanzuckerberg.com>.

## Configuration options

The following configuration options require special consideration in any multi-user or hosted environment:

`--disable-diffexp`: the differential expression computation can be resource intensive, in particular for large datasets. If many differential expression calculation requests are made in rapid sequence, it may cause the server CPU or memory resources to be exhausted, and impact the ability of other users to access data. This command line option will disable the differential expression feature, including the removal of the `Differential expression` button.

`--experimental-annotations`: this feature, which is disabled by default, may not be appropriate for hosted environments. It will write to the local file system, and in extreme cases could be used to abuse (or exceed) file system capacity on the hosting server.

`--experimental-annotations-file`: this specifies a single file for all end-user annotations, and is incompatible with hosted or multi-user use of cellxgene. Using it will cause loss of user annotation data (ie, the CSV file will be overwritten). If you wish to explore using the _experimental_ annotations feature in a multi-user environment, please refer to the [annotations documentation](annotations).

## Community software projects

There are a number of teams building tools or infrastructure to better utilize cellxgene in a multiple user environment. While we do not endorse any particular solution, you may find the following helpful.

- [Novartis Cellxgene Gateway](https://github.com/Novartis/cellxgene-gateway) - a multiple-user and multiple-dataset gateway for cellxgene.
- Interactive Enviroment in the [Galaxy Project](https://galaxyproject.org/) ([patch notes](https://docs.galaxyproject.org/en/release_19.05/releases/19.05_announce.html))

If you know of other solutions, drop us a note and we'll add to this list.

# Deploying cellxgene with Heroku

## Quickstart

Clicking on the following button will forward you to Heroku to begin the deployment process:

<a href="https://heroku.com/deploy?template=https://github.com/chanzuckerberg/cellxgene">
 <img src="https://www.herokucdn.com/deploy/button.svg" alt="Deploy">
</a>

If not already logged in to Heroku, there you will be prompted to log in or sign up for an account.

Once logged in you will be sent to the setup page. Here you can set some of the basic settings for the app:

### Default settings

- `App name`: the unique name for your deployment
- This will also serve as the default URL (e.g. https://cellxgene.herokapp.com/)
- `App owner`: Who will own this app. Either you personally or an organization/team
- `Region`: Location of the server where the app will be deployed (EU or US)

### Configuration

- `DATASET`: A _publicly_ accessible URL pointing to a .h5ad file to view
- This defaults to pbm3k.h5ad

After filling out the settings and pressing the `Deploy app` button Heroku will begin building your deployment. This process will take a few minutes, but once completed you will have a personal free hosted version of cellxgene!

## What is Heroku?

Heroku is a quick and easy way to host applications on the cloud.

A Heroku deployment of cellxgene means that the app is not running on your local machine. Instead, the app is installed, configured, and run on the Heroku servers (read: cloud).

On Heroku's servers, applications run on a [dyno](https://www.heroku.com/dynos) which are Heroku's implementation and abstraction of containers.

Heroku is one of many options available for hosting instances of cellxgene on the web.
Some other options include: Amazon Web Services, Google Cloud Platform, Digital Ocean, and Microsoft Azure.

## Why use Heroku to deploy cellxgene?

What Heroku enables is a quick, non-technical method of setting up a cellxgene instance. No command line knowledge needed. This also allows machines to access the instance via the internet, so sharing a visualized dataset is as simple as sharing a link.

Because cellxgene currently heavily relies on its Python backend for providing the viewer with the necessary data and tooling, it is currently not possible to host cellxgene as a static webpage.

This is a good option if you want to quickly deploy an instance of cellxgene to the web. Heroku deployments are free for small datasets up to around 250MBs in size. See below regarding larger datasets.

## When should I not deploy with Heroku?

- The default free dyno offered by Heroku is limited in memory to 512 MBs
  - The amount of memory needed for the dyno is roughly the same size as the h5ad file
  - Heroku offers tiered paid dynos. More can be found on the [Heroku pricing page](https://www.heroku.com/pricing)
  - Note that this can get _very_ expensive for larger datasets (\$25+ a month)
- On the free dyno, after 30 minutes of inactivity, Heroku will put your app into a hibernation mode. On the next access, Heroku will need time to boot the dyno back online.
- Having multiple simultaneous users requires more memory. This means that the free container size is easily overwhelmed by multiple users, even with small datasets; this can be addressed by purchasing a larger container size
- For this facilitated Heroku deployment to work, your dataset must be hosted on a publicly accessible URL
- By default, Heroku publically shares your instance to anyone with the URL.
  - There are many ways of securing your instance. One quick and simple way is by installing [wwwhisper](https://elements.heroku.com/addons/wwwhisper), a Heroku addon
