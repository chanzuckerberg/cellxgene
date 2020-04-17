FROM busybox

RUN echo Hello World from $(hostname) > index.html
CMD busybox httpd -f -p ${PORT}

