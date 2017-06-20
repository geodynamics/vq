#!/bin/bash

groupmod --gid=$HOST_GID virtualquake
usermod --uid=$HOST_UID virtualquake
chown -R virtualquake:virtualquake /home/virtualquake/
cd /home/virtualquake/external_volume
su virtualquake
