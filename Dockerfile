# FROM busybox
FROM python:3.6

RUN apt update -y && apt upgrade -y
RUN apt install -y make tree git build-essential
# RUN echo Hello World from $(hostname) > index.html
# CMD busybox httpd -f -p ${PORT}

ENV ENV=dev

RUN mkdir -p /app/cellxgene
WORKDIR /app/cellxgene

# TODO(el): copy everything for now... figure out later
COPY server server
RUN pip install -e .

CMD python server/eb/app.py
# CMD make build-client
# CMD cd server/eb && make build
