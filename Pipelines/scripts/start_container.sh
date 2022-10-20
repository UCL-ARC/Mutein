#!/bin/bash

#ensure unique tmp folders for each container instance
if [ -z "${TMPDIR+x}" ]; then
    #TMPDIR is not set, use a tmp under mutein working directory
    CONTAINER_DIR_BASENAME=${MUT_DATA}/tmp/singularity.${RANDOM}${RANDOM}
else 
    #TMPDIR is set, use fast local disk
    CONTAINER_DIR_BASENAME=${TMPDIR}/singularity.${RANDOM}${RANDOM}
fi

export SINGULARITY_CACHEDIR=${CONTAINER_DIR_BASENAME}/cachedir
export SINGULARITY_TMPDIR=${CONTAINER_DIR_BASENAME}/tmpdir
export SINGULARITY_LOCALCACHEDIR=${CONTAINER_DIR_BASENAME}/localcachedir
export SINGULARITY_PULLFOLDER=${CONTAINER_DIR_BASENAME}/pullfolder

mkdir -p ${SINGULARITY_CACHEDIR} \
         ${SINGULARITY_TMPDIR} \
         ${SINGULARITY_LOCALCACHEDIR} \
         ${SINGULARITY_PULLFOLDER}

#run a bootstrap script which ends with "bash" to give an interactive session?
#need to set PS1 it seems?! may need to set home somehow and not map cwd??
singularity shell --no-home --bind /lustre containers/mutein.sif

# one way to run a setup script in the container but still get
# an interactive session
##singularity shell containers/mutein.sif \
##  -c "bash script_singularity.sh && /bin/bash -rcfile singularity_bashrc"