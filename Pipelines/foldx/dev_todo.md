---------------------------------------------------------------------------
#### dev todo list
- Add a param for the user's home directory in the main bash scripts
- Make the number of repairs confugurable and use the repaired file with the specified number from thruputs
- Get python actions for continuous integration working
- Get a so called "empty" run that can work in CI (foldx won't work I assume???)
- Make the code and comments and variables up to rsdg standard
- FOR CI - add a final test to add asserts for changed data (waiting to be happy with data#1)
- If one of the array jobs fails the dependent batch can still go ahead. Do I just re-run the jobs that failed (they are numbered).
- the array number for the variant split is the combination of individual variants, so it should be calculated rather than passed in
---------------------------------------------------------------------------
#### discuss todo list 
- Anaytics - establish agreement or what the future action is on the number differences as per above
---------------------------------------------------------------------------
### In progress
- Make the batch dependency, script names and times (everything for qsub) configurable from file not hard coded
- complete the pipeline for the variants in the style adopted (batches 5,6,7)

### done
#### sprint 1 - 7/4/22
DEV
- Put the repairs in another folder
- Log the input paramaters to the current directory
- Implement a full dummy test data set to mimic an HPC batch
- Change the files so there is only 1 python file per bash script (__main__)
- make the test set 1tst for potential multiple ones
- Make a fully demonstrable show of the orders-of-mutations problem in its own test file
DISCUSS
