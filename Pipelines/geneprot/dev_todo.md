---------------------------------------------------------------------------
### Before a pull request
- Run the pep8 code formatter (https://pypi.org/project/black/)
```
python -m black Pipelines/geneprot/scripts/*.py
python -m black Pipelines/shared/lib/*.py
```
- Run the manual regression tests

#### dev todo list


- Make a lib higher up to be used by all pipelines
- A new structure is needed to be used by all pipelines, or an outer and inner file structure for the pipelines, ie a gene could have a folder with multiple pdbs, one of those pdbs may be hosen
- work out how to do the web scraping, the pdb just announced a new interface which would definitely be the best thing to do: http://search-beta.rcsb.org/#search-api

