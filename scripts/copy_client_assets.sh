#!/bin/bash

set -e 

# script to copy the client static assets to a location known to
# the server.  

if [ $# != 2 ] ; then
   echo "Usage $0 <client> <server>"
   exit 1
fi

CLIENT=$1
SERVER=$2

mkdir -p $SERVER/common/web/static/img
mkdir -p $SERVER/common/web/static/js
mkdir -p $SERVER/common/web/templates/
cp $CLIENT/index.html $SERVER/common/web/templates/
cp -r $CLIENT/static $SERVER/common/web/
cp $CLIENT/favicon.png $SERVER/common/web/static/img
cp $CLIENT/service-worker.js $SERVER/common/web/static/js/

