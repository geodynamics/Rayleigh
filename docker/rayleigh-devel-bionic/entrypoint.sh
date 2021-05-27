#!/bin/sh

groupadd -g $HOSTGID $HOSTUSER
useradd -u $HOSTUID -g $HOSTGID -d /work -s /bin/bash $HOSTUSER

echo "$HOSTUSER ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/nopasswd

su - $HOSTUSER
