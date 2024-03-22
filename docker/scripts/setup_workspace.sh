#!/bin/bash

#
# setup the pipeline workspace folders
# run from inside the container once running inside the TRE
# run from inside the docker subfolder
#

source settings

#create the raw and workspace+config subfolders of the main data directory
mkdir -p "${GUEST_DATA_FOLDER}/mutein_raw"
mkdir -p "${GUEST_DATA_FOLDER}/mutein_workspace/config"

#create the site-specific settings file from the template file
cp "/home/${CONTAINER_USER}/${REPO_FOLDER}/Mutein/mutein/config/mutein_settings_tre" "${GUEST_DATA_FOLDER}/mutein_workspace/config/mutein_settings"

#fill in the site-specific values of the settings file
cd "${GUEST_DATA_FOLDER}/mutein_workspace/config"
export NEWVAL="/home/${CONTAINER_USER}/${REPO_FOLDER}/Mutein" ; sed -i "s|MUT_DIR_PLACEHOLDER|${NEWVAL}|"  mutein_settings
export NEWVAL="${GUEST_DATA_FOLDER}/mutein_workspace" ;         sed -i "s|MUT_DATA_PLACEHOLDER|${NEWVAL}|" mutein_settings
export NEWVAL="${GUEST_DATA_FOLDER}/mutein_raw" ;               sed -i "s|MUT_RAW_PLACEHOLDER|${NEWVAL}|"  mutein_settings


