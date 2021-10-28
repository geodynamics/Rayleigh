#!/bin/sh

export HOME=/work
cd "$HOME"

if [ -n "$HOSTGID" -a -n "$HOSTUSER" -a -n "$HOSTUID" ]
then
	# only create user if all options are provided

	# This only creates the group if a group with that GID does not already exist.
	getent group $HOSTGID > /dev/null || groupadd -g $HOSTGID $HOSTUSER
	useradd -u $HOSTUID -g $HOSTGID -d /work -s /bin/bash $HOSTUSER

	echo "$HOSTUSER ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/nopasswd

	sudo -u $HOSTUSER /bin/bash --rcfile /etc/bash.bashrc -i
else
        if [ "${NOUIDWARN}" != 1 ]
        then
                echo
                echo "You are running this container without providing the HOSTUSER, HOSTGID, and HOSTUID"
                echo "environment variables. This means the container will use the root user and it might create files"
                echo "owned by root in your home directory if you run this on Docker with Linux as a host system."
                echo 'Consider running this using the "docker-devel" script from the Rayleigh main directory.'
                echo
        fi
	/bin/bash --rcfile /etc/bash.bashrc -i
fi
