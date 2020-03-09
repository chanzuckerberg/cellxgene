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

1. Create an S3 bucket

   Upload your matrix files to this bucket
   
2. Create an elastic beanstalk application.  For example:

    ```
    EB_APP=cellxgene-app
    eb init -p python-3.6 $EB_APP
    ```

3. Create the artifact.zip file for the application

   ```
   make build
   ```
   
4. Create an environment

    ```
    # name of the environment
    EB_ENV=cellxgene-env
    # type of ec2 instance to run the cellxgene server.     
    EB_INSTANCE=m5.large 
    CXG_DATAROOT=<location to your S3 bucket>
   
    eb create $EB_ENV --instance-type $EB_INSTANCE --envvars CXG_DATAROOT=$CXG_DATAROOT
    ```

5. Give the elastic beanstalk environment access to the S3 bucket.

    This link may provide some useful information:
    https://aws.amazon.com/premiumsupport/knowledge-center/elastic-beanstalk-s3-bucket-instance/
    
6. Deploy the application

    ```
    eb deploy $EB_ENV 
    ```
   
7. Open the application in a browser

   ```
   eb open
   ```
