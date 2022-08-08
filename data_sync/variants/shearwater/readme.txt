Inigo has now confirmed we can share these scripts. The license is GPL-3. 

There are a few versions of these scripts because they have been adapted simultaneously by different research groups.  I am attaching the ones we have been sent by Inigo Martincorena’s group.

The particularly complicated part is the filtering. How this is done depends on the data set, so it may not be appropriate to just copy the scripts that have been used before. 

I’ve attached the three original scripts I was sent in 2017:
wrapper_shearwaterML_multiplebams_orig.R   - this runs the variant calling on some bait set regions. This script contains a small bug (see below). 
shearwater_pipe_WithPredefinedBAMList.R   - this is a script that loops through a bait set bed file and launches wrapper_shearwaterML_multiplebams_orig.R  on the Sanger compute farm. 
annotate_mutations.R  - although it is called annotate, it is more of a filtering script.  

I was also sent a couple of updated scripts in 2020. There was a minor bug in the wrapper script which is fixed here and there may be other changes. And the annotate_mutation script is updated for the project they were working on at the time:
annotate_mutations_vAL.R
wrapper_shearwaterML_multiplebams_v120719.R

