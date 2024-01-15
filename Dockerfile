# Use Ubuntu as the base image
FROM ubuntu:latest AS base

# Install build-essential, cmake, and OpenBLAS
RUN apt-get update && \
    apt-get install -y build-essential cmake libopenblas-dev

# Set environment variables
ENV OpenBLAS_HOME=/usr/include/x86_64-linux-gnu
ENV OpenBLAS=/usr/lib/x86_64-linux-gnu
