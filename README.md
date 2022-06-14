# Mutein
#### Overview
RSDG project for [Ben Hall](https://iris.ucl.ac.uk/iris/browse/profile?upi=BHALL50), original tracker page [here](https://github.com/UCL-RITS/research-software-opportunities/issues/549). This repository is intended to stay private within ARC and the Hall lab and functions to track internal issues within the project. Any future public release of the source code, such as for a publication, should copy the code into a new, separate repository to avoid publishing the internal discussions. This repo acts as an onboarding resource for any future RSDG developers needing to work on the project. We use the github [issues](https://github.com/UCL/Mutein/issues) and [kanban board](https://github.com/UCL/Mutein/projects/1) as the main project management and logging tools.

This is a research project to reanalyse multiple existing [datasets](https://github.com/UCL/Mutein/issues?q=is%3Aissue+label%3Adataset) from publications, aiming to understand the effects of mutations on protein structure, and through that to ultimately understand the proteins' functions and interactions. The datasets are DNA sequence data, including target capture of specific gene sets, whole exome capture and whole genome capture, from studies looking at diseased and normal tissue from various tissue types.

The project will essentially scale up and improve upon what has already been done in individual papers on single datasets by Ben and/or other lab members or collaborators (such as in [Fowler et al 2021](https://github.com/UCL/Mutein/issues/54), and produce a well documented, maintainable pipeline that Hall lab members can reuse or build on in the future.

The data are being downloaded and processed on UCL's [Myriad](https://www.rc.ucl.ac.uk/docs/Clusters/Myriad/) HPC.

#### Data Analysis Pipelines
The data analysis process involves cloning this repository, installing the required dependencies, then running the appropriate scripts. Due to various sources of [unreliability](https://github.com/UCL/Mutein/issues/95) it probably not feasible to run absolutely the entire pipeline from start to finish in a fully automated fashion. Never the less it is a key aim to have all the analysis steps scripted and automated as far as possible, and to support the use of command line scripting to accomplish every step, due to the suitability of scripting for self documentation, automation, maintainability and scaleability, and to enhance the shareability and reproducibility of the entire analysis

