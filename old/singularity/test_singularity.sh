#!/bin/bash

#
# run a pipeline 
# that tests the singularity features of yamlmake
#
# run on myriad using:
# cd [project working directory]
# source config/mutein_settings
# ${MUT_DIR}/Pipelines/singularity/test_singularity.sh
#

yamlmake --yaml ${MUT_DIR}/Pipelines/yamlmake/testing/test_singularity.yml --dryrun
