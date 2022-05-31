# MuteinGUI

## disclaimer
I have never used tkinter before so this is a very rough first attempt.

## To pip install from vscode, find your python, e.g.
```
C:\Users\Rache\AppData\Local\Programs\Python\Python310\python.exe -m pip install seaborn
```

# Running the MuteinGUI

The MuteinGI is in the Mutein github in the MuteinGUI directory. It is compeltely standalone from the rest of Mutein.

It is a pure client, accessing the servers using ssh calls. Via ssh calls it submits qsub or python scripts.

All actions done by the GUI can also be done in the standard bioinformatics way of ssh-ing onto the server and running the approriate scripts. Full documentation is given for this here (in this currently doesn't exist document).

Once installed you can run MutinGUI by either
- Navigating to it in an idea and running the file muteingui.py
- Running from a commandprompt 
->python muteingui.py
-- (or full paths for both as necessary)

The MuteinGUI currently assumes that Mutein has been installed as per these instructions: (link when the doc is updated)

## TABS 
### Config
- Input your username and password
- leave the server as myriad unless you are using a different one.
- Press regenerate paths, you will see the server paths
- Press Verify connection to check all is well

### Monitor
This is the main page to use to monitor th state of jobs on queues. The buttons are importantly coloured - think twice if the button is not green as an action is taken - deleting log files ar cancelling jobs or re-running them.
- Refresh QStat, just like running qstat on the box.
- Clean (orange) will delete all succesfully completed jobs (no errors, a complete tag) and flag error files.
- Delete all log files (red) only do this when your between jobs and nothing is running or being investigated.
- Rerun errors (red) although not too dramatic, it will resubmit any of the jobs with error files.
- View Log, thuis will use the entry box to match any log files with that string and show them.
-  Cancel Jobs (red) just one is fine but be very careful if you put an asterisk in as it will cancel all jobs.

### Pipeline
This is the pipeline page from where you can manage the pipeline submissions. You can look at the pipeline from a dataset, gene or pdb perspective. You may only be looking at a pdb in isolation in which casethere is nothing in the dataset or gene box, or you may be focusing on a pdb within a dataset, in which case those boxes will be filled. Likewise a gene.
- Submit pdb prepare
-- For a whole dataset it will (from vcf file ultimately but...) establish the genes in scope and their variants, and  then for each gene...
-- for a gene it will use swiss model and unitpor and alphafold to find all the structures and establish the coverage sections from gene to pdb
- Submit repair
-- each structure is run through the foldx energy repair program 10 times to optimise the energy of the structure
- Submit splits prepare
-- The residues and variants are chunked up into tasks so that the intensive calculations can be parallelised. This task runs in pythong as it is quick and avoids waiting for the q submission.
- Submit tasks (red)
-- This is in red because be warned that there could be many tasks and you could be waiting days/weeks for the results. Kicking off an entire dataset for example could be many thousands of 8 thread tasks. Just 1 pdb (unless a very large em structure) is likely to be reasoably quick).
- Submit missing task (red)
-- Only press this if you have brought down the currently running tasks. This will recognise those tasks that have not completed and resubmit just those. Mostly it is in case of a server failure/user error, it shouldn't be necessary normally. 
- Submit aggregation
-- At the gene level a summary of all the structures are put in a single dataframe. At the pdb level all the parallelised tasks are summarised into 1 dataframe. In fact 3 dataframes as the background ddg is in 1, and then the position scan and build model ddg for the variants in 2 more.

### Results
Not currently complete, this is just a beginning, it shows the 3 dataframes in a text box.







    






