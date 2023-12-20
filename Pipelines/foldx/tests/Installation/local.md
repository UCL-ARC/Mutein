# Test an installation on your local machine
- (Also provides documentation)
----
- Add yourself to the environments.py in libs

### The data has 3 levels
- 1 dataset has many genes
- 1 gene has many protein structures (pdbs)
- 1 pdb is split into many tasks for parallelisation

- The folder


### There are 3 tests scripts to run to verify local installation
- pdb, this tests the main functionality
- gene, this additionally tests finding protein structures and aggregating them all to a final file
- dataset, this additionally handles many genes

- Run the small test script in this directory install_tst.py
