#!/bin/bash

export PATH=${PATH}:/opt/singularity/bin
sudo singularity build --allow-setuid mutein.sif mutein.def
