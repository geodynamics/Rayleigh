#!/bin/bash

# This script generates a docker image from the latest Rayleigh development version.
# It requires a docker installation on the local machine, and the ability to
# communicate with the docker daemon without root user privileges (see the docker
# webpage for an explanation).
# Note: This container is build from the developer version on Github, it does not use
# the local Rayleigh folder. Therefore local changes are not included in the container.

docker build --no-cache -t geodynamics/rayleigh:tacc-latest . 2>&1 | tee build.log
