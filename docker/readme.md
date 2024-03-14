# TRE docker container

The folder contains scripts used to setup a docker container to run the Mutein pipeline in a TRE environment. The container could be created on any machine capable of building docker images, but is being developed on a Vagrant managed Ubuntu VM running under VirtualBox on a Windows 10 machine.

## Summary
- The container needs to be setup to contain all the dependencies to run the pipeline from inside the TRE
- A cutdown subset of the data and associated metadata will be downloaded to this local VM and uploaded to the TRE separately using an scp push

## Instructions
- [install_docker.sh](scripts/install_docker.sh) used to install docker on the Ubuntu VM
