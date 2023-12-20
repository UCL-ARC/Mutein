#!/bin/bash

#
# this script is to install and build singularity as root on a local vm you have admin rights for
# this is to allow your local version of singularity to run images without first converting them to sandbox
# to address to issue described here: https://github.com/apptainer/singularity/issues/6065
# otherwise the built container must be (labouriously) converted to a sandbax locally for testing
# or else uploaded to myriad for testing
#

# see here for details: https://docs.sylabs.io/guides/3.5/admin-guide/installation.html
# but those instructions didn't work
# so I followed the INSTALL.md instructions I found in the downloaded tarball

#install dependencies
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    uuid-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup-bin

sudo apt-get update && \
  sudo apt-get install -y build-essential \
  libssl-dev uuid-dev libseccomp-dev \
  pkg-config squashfs-tools cryptsetup

#also needs go-lang >= 1.13
#installed golang-1.13 using golang metapackage on ubuntu 20.04 LTS
sudo apt install -y golang

#clone the repo
mkdir -p ${HOME}/go/src/github.com/sylabs
cd ${HOME}/go/src/github.com/sylabs
git clone https://github.com/sylabs/singularity.git
cd singularity
git checkout v3.5.3

#build singularity as root
./mconfig --prefix=/opt/singularity
cd ./builddir
sudo make #note: trying to enable running images without sandboxing by compiling as root
sudo make install

