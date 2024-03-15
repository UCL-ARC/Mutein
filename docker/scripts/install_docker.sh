#!/bin/bash

#
# script was run on a Vagrant managed Ubuntu 22.04 VM
# to allow the TRE container to be built on the docker build VM
# prior to uploading it to the TRE
#

# see also
# https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-22-04

sudo apt update
sudo apt install -y apt-transport-https ca-certificates curl software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt update
apt-cache policy docker-ce
sudo apt install -y docker-ce
sudo systemctl status docker

#allow sudoless access to docker commands
sudo adduser ${USER} docker
#now log out and in again
