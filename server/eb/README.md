# AWS Elastic Beanstalk

This directory contains script to aid in creating and deploying cellxgene on
an AWS Elastic Beanstalk instance.

This will result in a variant of cellxgene, running on AWS EC2 instances, serving data from S3.
All datasets must be in the new CXG (tiledb) format - see the converter script cxgtool.py
in server/converters - and located in a single S3 prefix, which is accessible to the instance.
In the current incarnation, no access control or authentication support is available
(outside of anything you configure yourself), so this is most appropriate for public datasets.

This is early development work, and will change significantly in the near future.
We would love feedback on it, but please assume it will change.

## Prerequisites

1. Some familiarity with AWS EB, S3, and IAM are needed.

2. Install the awsebcli.
Instruction are here:
https://docs.aws.amazon.com/elasticbeanstalk/latest/dg/eb-cli3-install.html

3. In the top level directory, run ```make build-client``` to create the client static assets.

## Steps

These steps are meant to serve as an example.
There are many more options to these commands that may be important or necessary for your environment.

### 1. Make your matrix files available to the EB servers.

The following choices are known to work.

* S3 Bucket.
* POSIX filesystem (such as Lustre)
* Lustre filesystem backed by S3

S3 is convenient and the relatively inexpensive option.
Lustre is higher performance, but more expensive, and slightly more complex to setup and manage.
AWS supports a feature to back the Lustre filesystem with S3, which give an easy to manage and high
performance option.

Once the storage is in place, the next step is to copy your matrix files to that location.
Currently cellxgene supports a flat file organization.   Each matrix file is located from
the same s3 prefix or filesystem directory.  This location is specified in the configuration as the dataroot.

### 2. Create an elastic beanstalk application.  For example:

```
EB_APP=cellxgene-app
eb init -p python-3.6 $EB_APP
```

### 3. Configuring cellxgene

All the cellxgene configuration options can be set from a configuration file.
This file can be generated like this:

```cellxgene launch --dump-default-config > myconfig.yaml```

The config file may then be customized before the app is deployed.

There are two ways to set the config file location, evaluated in this order:

First, if your config file is named "config.yaml" and exists in `customize/config.yaml`,
then it will be bundled with the application zip file and installed along
side the app on the EB servers.

Second, a potentially more flexible approach is to place your config file in a location accessible to the EB
servers, such as in S3.  For example:  s3://my-bucket/my-datasets/config.yaml.
Set the CXG_CONFIG_FILE environment variable to specify this location.

Another option is to set the CXG_DATAROOT environment variable.  The dataroot
is the location where the matrix files are located.
This environment variable will override the dataroot in the config file (if specified).

- Note:  Certain features, such as user annotations, are automatically disabled by the EB app,
and cannot be enabled using configuration.  They may be enabled manually by modifying app.py, however
this is not supported or recommended at this time.

### 4. Customization

The deployment can be customized in several ways, by adding files to a directory called
`customize` which is placed in this directory.

#### config file

This was described in the previous section.

#### static files

The cellxgene server can serve additional static webpages that will be associated with the app.
These include the about_legal_tos (terms of service), and about_legal_privacy, for example.
To use this feature, do the following:

* In this directory, create a sub directory called "customize/deploy/".
* Copy the files you want to serve into this directory
* modify your configuration file to set the location to these file:  /static/cellxgene/deploy/<filename>

Example:  you want to include an "about_legal_tos" and "about_legal_privacy" page to cellxgene.
Assume files called "tos.html" and "privacy.html" exist.

```
$ mkdir static
$ cp <source_dir>/tos.html customize/deploy/tos.html
$ cp <source_dir>/privacy.html customize/deploy/privacy.html

# edit config.yaml
$ grep "/static/cellxgene/deploy" config.yaml
about_legal_tos: /static/cellxgene/deploy/tos.html
about_legal_privacy: /static/cellxgene/deploy/privacy.html
```

#### Inline javascript scripts

Additional scripts can be added using the server/inline_scripts config parameters.
To include these scripts in the deployment, use the following steps:

* In this directory, create a sub directory called "customize/inline_scripts".
* Copy the script files into this directory
* modify your configuration file to set the location to these file (leaving off customize/inline_scripts)

For example, to add an inline script called "myscript.js":

```
$ mkdir scripts
$ cp <source_dir>/myscript.js customize/inline_scripts/myscript.js
# edit the config.yaml
$ grep inline_scripts config.yaml
  inline_scripts : [ myscript.js ]
```

#### Plugins

Optionally, you can add plugins to the server python code. To include a plugin in the deployment use the following steps:

```
$ mkdir plugins
$ cp <source_dir>/<my_plugin>.py customize/plugins/<my_plugin>.py
```

#### ebextensions

Any additional config files intended for the `.ebextensions` directory of the artifact can be added
to the `customize/ebextensions` directory.  Any file found here will be copied over.

#### requirements.txt

A custom requirements.txt can be supplied in customize/requirements.txt.
This file must fully specify the versions of all the python modules used by the server in the deployment.
This is useful to ensure that the dependencies do not change from one deployment to the next.
Therefore the custom/requirements.txt must all have exact versions specified (e.g. anndata==0.7.1).

This file can be generated the first time using a process like this:
```
#  assume you are running in this directory
$  virtualenv temp
$  source temp/bin/activate
$  pip install -r ../requirements.txt
$  mkdir -p customize
$  pip freeze > customize/requirements.txt
$  deactivate
$  rm -rf temp/
```

Keep the customize/requirememts.txt file, and reuse it for each deployment.
If a future cellxgene version updates its requirements by modifying a module version
or adding a new dependency, then the `make build` process will detect any
incompatibilities and raise an error.

### 5. Create the artifact.zip file for the application

```
$ make build
```

### 6. Flask secret key

The application requires as secret key to be provided to flask, the web framework used by cellxgene.
There are three ways to provide the secret key:

- In the configuration file: update the server/flask_secret_key attribute.
- An environment variable: `CXG_SECRET_KEY`
- Managed by the AWS Secret Manager

If using the AWS Secret Manager, then the secret name is passed as an environment variable:  CXG_AWS_SECRET_NAME.
The secret must contain a key with the name "flask_secret_key".
The region name for the AWS Secret Manager must be specified (e.g. us-east-1).
The most straightforward way is to specified it with the CXG_AWS_SECRET_REGION_NAME environment variable.
If this environment variable is not defined, then the app attempts to determine the region from the
dataroot (if in s3), or the config file location (if in s3).

### 7. Create an environment

```
# name of the environment
$ EB_ENV=cellxgene-env

# type of ec2 instance to run the cellxgene server (for example)
$ EB_INSTANCE=m5.large

# One or both of the following environment variables needs to be set
$ CXG_DATAROOT=<location to your S3 bucket>
$ CXG_CONFIG_FILE=<location to your config file>

# Potentially also set envvars for the secret key.

$ eb create $EB_ENV --instance-type $EB_INSTANCE \
   --envvars CXG_DATAROOT=$CXG_DATAROOT,CXG_CONFIG_FILE=$CXG_CONFIG_FILE
```

### 8. Give the elastic beanstalk environment access to the dataroot.

If using S3, this link may provide some useful information:
https://aws.amazon.com/premiumsupport/knowledge-center/elastic-beanstalk-s3-bucket-instance/
If using Lustre, then this link may provide a place to start:
https://aws.amazon.com/fsx/lustre/

### 9. Deploy the application

```
$ eb deploy $EB_ENV
```

### 10. Open the application in a browser

```
$ eb open $EB_ENV
```
