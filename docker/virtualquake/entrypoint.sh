#!/bin/bash

if [ -n "${HOST_UID}" ]
  then
    groupmod --gid=$HOST_GID virtualquake
    usermod --uid=$HOST_UID virtualquake
fi
chown -R virtualquake:virtualquake /home/virtualquake/
cd /home/virtualquake/external_volume
su virtualquake
