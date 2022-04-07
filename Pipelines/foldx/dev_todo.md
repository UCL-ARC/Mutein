---------------------------------------------------------------------------
#### dev todo list
- I am only being very simple in job splits - can foldx multithread?
- Make the code and comments and variables up to rsdg standard
- FOR CI - add a final test to add asserts for changed data (waiting to be happy with data#1)
- Make the CI empty the data completely first
- If one of the array jobs fails the dependent batch can still go ahead. Do I just re-run the jobs that failed (they are numbered).
- the array number for the variant split is the combination of individual variants, so it should be calculated rather than passed in
- unsatisfactory use of both split and array jobs needs thinking about
- use argparse (https://docs.python.org/3/library/argparse.html#)
- Make foldx a "runner" so it is easier to understand what has gone into it.
---------------------------------------------------------------------------
#### discuss todo list 
- Anaytics - establish agreement or what the future action is on the number differences as per above
---------------------------------------------------------------------------
#### sprint 1 - 7/4/22
DEV


DISCUSS
- 

DEV-ED
- Put the repairs in another folder
- Log the input paramaters to the current directory
- Implement a full dummy test data set to mimic an HPC batch
- Change the files so there is only 1 python file per bash script (__main__)
- make the test set 1tst for potential multiple ones
- Make a fully demonstrable show of the orders-of-mutations problem in its own test file
- Make the batch dependency, script names and times (everything for qsub) configurable from file not hard coded
- Make the number of repairs confugurable and use the repaired file with the specified number from thruputs
- Make readme.md the main file "information that the next generation of devs"
- enable override of the pieline params from the batch, eg combos array num changes based on pdb size it is not an arbitrary 
performance decision)
- complete the pipeline for the variants in the style adopted (batches 5,6,7)
- Get rid of magic number in get_make_paths
- Make args a class to handle args (use argsparse too - later)
- Add a param for the user's home directory in the main bash scripts
- Get python actions for continuous integration working
- Get a so called "empty" run that can work in CI (foldx won't work I assume???)

DISCUSSED-ED
- We need to change the analysis for pos scan due to bug
- the 10x repair is more stable with the ddg for the histidine, have results to show, so yes but there is also a problem
- Another bug in the copying across of the actual results to results file is why the build script is a bit weird
---------------------------------------------------------------------------
