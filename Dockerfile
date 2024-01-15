# Use Ubuntu as the base image
FROM ubuntu:latest AS base

# Install build-essential, cmake, and OpenBLAS
RUN apt-get update && \
    apt-get install -y build-essential cmake libopenblas-dev
