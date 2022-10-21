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

# singularity exec --no-home --bind /lustre containers/mutein.sif \
#     "/usr/bin/bash -l ${MUT_DIR}/Pipelines/scripts/setup_container.sh"

#--contain removes fuse device from /dev
#need to prevent /tmp:/tmp mapping without losing /dev/fuse

singularity exec \
    --no-home --contain --writable-tmpfs \
    --bind /lustre \
    containers/mutein.sif \
    ${MUT_DIR}/Pipelines/scripts/setup_container.sh \
    ${MUT_DATA}/cipher