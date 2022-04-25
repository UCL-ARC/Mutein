---------------------------------------------------------------------------
### Before a pull request
- Run the pep8 code formatter (https://pypi.org/project/black/)
```
python -m black Pipelines/geneprot/scripts/*.py
python -m black Pipelines/shared/lib/*.py
```
- Run the manual regression tests

#### dev todo list
- The override of the time or array for a batch needs to make sense
- need to make sure the covid batch still works on hpc
- need to scope out how long things could take
- I might need to change the scripts now to be HPC-splittable on a per-gene basis?


-------------------------------------------------------------------
Sprint #3
TODO
- Does any of it work on hpc now?
- Make it so I can choose to use mouse or human or other and all that entails including the alphafold structures
- Change the analysis for the foldx pipeline to NOT use combinations
- Do a full test on notch1, human and mouse
- The tests no longer work 
DONE

-------------------------------------------------------------------
Sprint #2
DONE
- Make a lib higher up to be used by all pipelines
- A new structure is needed to be used by all pipelines, or an outer and inner file structure for the pipelines, ie a gene could have a folder with multiple pdbs, one of those pdbs may be hosen
- work out how to do the web scraping, the pdb just announced a new interface which would definitely be the best thing to do: http://search-beta.rcsb.org/#search-api - I have used bioservices
- comments before pull request
- pep8 before pull request
- I have shearwater hardcoded at the moment

