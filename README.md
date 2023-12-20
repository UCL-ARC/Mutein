# Mutein
#### Overview
RSDG project for [Ben Hall](https://iris.ucl.ac.uk/iris/browse/profile?upi=BHALL50), original tracker page [here](https://github.com/UCL-RITS/research-software-opportunities/issues/549). This repository is intended to stay private within ARC and the Hall lab and functions to track internal issues within the project. Any future public release of the source code, such as for a publication, should copy the code into a new, separate repository to avoid publishing the internal discussions. This repo acts as an onboarding resource for any future developers needing to work on the project.

This is a research project to reanalyse multiple existing [datasets](https://github.com/UCL/Mutein/issues?q=is%3Aissue+label%3Adataset) from publications, aiming to understand the effects of mutations on protein structure, and through that to ultimately understand the proteins' functions and interactions. The datasets are DNA sequence data, including target capture of specific gene sets, whole exome capture and whole genome capture, from studies looking at diseased and normal tissue from various tissue types.

The project will essentially scale up and improve upon what has already been done in individual papers on single datasets by Ben and/or other lab members or collaborators (such as in [Fowler et al 2021](https://github.com/UCL/Mutein/issues/54)), and produce a well documented, maintainable pipeline that Hall lab members can reuse or build on in the future.

#### Preparation
In order to use an HPC to accelerate the pipelines you will need to be setup with an account on the appropriate system(s). In the UCL context this currently means [Myriad](https://www.rc.ucl.ac.uk/docs/Clusters/Myriad/) for non-sensitive, non-GDPR data, and [CS](https://hpc.cs.ucl.ac.uk/) ([contact](https://hpc.cs.ucl.ac.uk/contact-us/)) for sensitive data.

See [working_on_cs.md](variant_calling/working_on_cs.md) for how to get started working on the CS HPC system with Mutein.

