import subprocess
import pytest
import os
import datetime
from datetime import timezone
from collections import namedtuple

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

# Starting from an empty data dir, running the failing pipeline repeatedly
# should keep stopping at the failure
@pytest.mark.parametrize("num_reps", [1, 2])
def test_idempotence_from_empty(num_reps, tmpdir):
    assert_count_items_dir(tmpdir, 0)
    for _repeat in range(num_reps):
        # run the workflow more than once to check that if output files were created by a
        # failed command, on subsequent attempts to run the workflow those files are not used and the
        # workflow continues to fail
        cp = run_yamlmake('basictest.yml', tmpdir)
        # should detect failure every time
        assert cp.returncode == 2
        assert os.path.isfile(tmpdir / 'good.txt')
        # file created by the failed command should still be present (might be needed for debugging) but
        # will be marked as bad by setting the special mtime value
        assert os.path.isfile(tmpdir / 'bad.txt')
        assert STALE_DATETIME == os.path.getmtime(tmpdir / 'bad.txt')
        # no unexpected files have crept in
        assert os.path.isdir(tmpdir / 'yamlmake_logs')
        assert_count_items_dir(tmpdir, 3)
        assert not os.path.exists(tmpdir / 'verybad.txt')

def test_run_until(tmpdir):
    assert_count_items_dir(tmpdir, 0)
    cp = run_yamlmake('basictest.yml', tmpdir, ['--run-until', 'jeremy'])
    # we instructed it to stop before the task that fails, therefore should return success
    assert cp.returncode == 0
    assert os.path.isfile(tmpdir / 'good.txt')
    assert os.path.isdir(tmpdir / 'yamlmake_logs')
    assert not os.path.isfile(tmpdir / 'bad.txt')
    assert not os.path.exists(tmpdir / 'verybad.txt')
    assert_count_items_dir(tmpdir, 2)
