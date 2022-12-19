import os
import subprocess
import shutil
import datetime
from datetime import timezone

# special mtime to indicate that the command that generated a file returned an error state
STALE_DATETIME = datetime.datetime(1980, 1, 1, tzinfo=timezone.utc).timestamp()

def assert_count_items_dir(tmpdir, expected_items_count):
    assert os.path.isdir(tmpdir)
    assert expected_items_count == len(os.listdir(tmpdir))

def find_test_yamlfile(base_yamlmake_file):
    return os.path.join(os.path.dirname(__file__), base_yamlmake_file)

# Run yamlmake as a separate process instead of as a python module, to
# make sure initialisation is properly performed
def run_yamlmake(base_yamlmake_file, working_dir, extra_args=None):
    if extra_args is None:
        extra_args = []
    yamlfilepath = find_test_yamlfile(base_yamlmake_file)
    newenv = os.environ.copy()
    newenv['TEST_DATA_OUTPUT'] = str(working_dir)
    return subprocess.run(['yamlmake', '--yaml', yamlfilepath] + extra_args, cwd=working_dir, env=newenv)

# Set up a certain pre-existing pipeline state so that a test can be run
def setup_existing_files(output_dir, setup_base_dir):
    # only work from an empty starting point
    assert_count_items_dir(output_dir, 0)
    setup_dir =  os.path.join(os.path.dirname(__file__), 'pre_existing_data', setup_base_dir)
    shutil.copytree(setup_dir, output_dir, dirs_exist_ok=True)
