# Creating a Mutein Gene Analysis (FoldX) Installation

#### First you need to make sure you have been added to the FoldX group.
- To do this you need to download an academic license from https://foldxsuite.crg.eu/licensing-and-services#academic-license
- Then you need to email rc-support@ucl.ac.uk and ask to be put in the FoldX groupon myriad, gicing proof of academic license.

When that has been done (a day or two) you can install Mutein Gene Analysis.

#### Steps for installation
- log on to the server, e.g.:
-- ssh -J ucbtlcr@socrates.ucl.ac.uk ucbtlcr@myriad.rc.ucl.ac.uk
------------------------
####  If you need to do so, generate an ssh key
- mkdir .ssh
- cd .ssh
- ssh-keygen
- (press enter all 3 times without typing anything)
- cat id_rsa.pub
- You can now copy this file, go to github website
- In your Settings/SSH and GPG keys, create one and paste it in
------------------------
#### Now you can clone the Mutein repo
- From the home root directory, cd ~ if you need to
- git clone git@github.com:UCL/Mutein.git
- cd into Mutein
- Checkout the correct branch if necessary, e.g.
-- git checkout rachel/foldxv1
------------------------
#### Now prepare the data input for the pipeline
- cd ~
- mkdir MuteinData
###### There is sample data in Mutein/data_sync/SmallDemos/
- cp -r ~/Mutein/data_sync/SmallDemos/dataset_cutdown/ ~/MuteinData/
- cp -r ~/Mutein/data_sync/SmallDemos/dataset_notch/ ~/MuteinData/
- cp -r ~/Mutein/data_sync/SmallDemos/pdb_1pb5/ ~/MuteinData/
------------------------
#### You are now ready to run a job!!!
###### The smallest job is the 1pb5 pdb which is a good test
- Open the GUI
- On the config page 
-- enter your username, password and server
-- Check the connections
- On the Pipeline page
-- delete the dataset and gene, and enter 1pb5 in the pdb box
-- Click on "submit repair"
------------------------------
#### Now monitor your batch
- On the Monitor tab
-- Press Refresh QStat and watch the progress


