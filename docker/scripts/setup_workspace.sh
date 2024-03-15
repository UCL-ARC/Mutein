#!/bin/bash

#
# script to run once the container is running inside the tre
# and can mount the external workspace folder
#

# todo: convert this dockerfile script into a normal bash script
# WORKDIR $data_path
# RUN \
#     mkdir -p $data_path/mutein_raw && \
#     mkdir -p $data_path/mutein_workspace/config && \
#     mkdir -p $data_path/somefoldername/subfolder
# RUN cp /home/$username/$repo_path/Mutein/mutein/config/mutein_settings_tre mutein_workspace/config/mutein_settings
# WORKDIR $data_path/mutein_workspace/config
# #RUN export NEWVAL=/home/$username/$repo_path/Mutein ; sed -i "s|MUT_DIR_PLACEHOLDER|$NEWVAL|" mutein_settings
