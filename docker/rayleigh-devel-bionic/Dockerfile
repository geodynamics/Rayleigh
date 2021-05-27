FROM geodynamics/rayleigh-buildenv-bionic

RUN apt update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt install --yes \
    sudo

ADD /entrypoint.sh /entrypoint.sh

ENTRYPOINT /entrypoint.sh
