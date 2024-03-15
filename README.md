# Mutein
ARC Collaborations project for [Ben Hall](https://iris.ucl.ac.uk/iris/browse/profile?upi=BHALL50), original tracker page [here](https://github.com/UCL-RITS/research-software-opportunities/issues/549). This repository is intended to stay private within ARC and the Hall lab and functions to track internal issues within the project. Any future public release of the source code, such as for a publication, should copy the code into a new, separate repository to avoid publishing the internal discussions. This repo acts as an onboarding resource for any future developers needing to work on the project.

This is a research project to reanalyse multiple existing [datasets](docs/dataset_papers.md) from publications, aiming to understand the effects of mutations on protein structure, and through that to ultimately understand the proteins' functions and interactions. The datasets are DNA sequence data, including target capture of specific gene sets, whole exome capture and whole genome capture, from studies looking at diseased and normal tissue from various tissue types.

The project will essentially scale up and improve upon what has already been done in individual papers on single datasets by Ben and/or other lab members or collaborators, and produce a well documented, maintainable pipeline that Hall lab members can reuse or build on in the future.

The variant calling pipeline in this repository is intended to be used in conjunction with the [MuteinFold](https://github.com/UCL-ARC/MuteinFold) pipeline. Variant calling identifies which somatic mutations are present in the samples, then MuteinFold can be used to estimate the effect of those mutations on protein folding. For documentation on MuteinFold see the [MuteinFold wiki](https://github.com/UCL-ARC/MuteinFold-Issues/wiki). These repositories are currently private so you will need to request access if you don't already have it.

For documentation of how to install and use the Mutein's variant calling pipeline see the [getting started](docs/getting_started.md) page.

The [docker](docker/readme.md) folder contains scripts used to setup a docker container to run the pipeline in a TRE environment.