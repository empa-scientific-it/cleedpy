FROM alpine:latest AS base

RUN apk update && \
    apk add --no-cache build-base cmake openblas-dev python3
