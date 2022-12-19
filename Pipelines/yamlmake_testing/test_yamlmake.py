from testing_shared import *
import pytest
import os

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
        assert os.path.isfile(tmpdir / 'good0.txt')
        assert os.path.isfile(tmpdir / 'good1.txt')
        # file created by the failed command should still be present (might be needed for debugging) but
        # will be marked as bad by setting the special mtime value
        assert os.path.isfile(tmpdir / 'bad.txt')
        assert STALE_DATETIME == os.path.getmtime(tmpdir / 'bad.txt')
        # no unexpected files have crept in
        assert os.path.isdir(tmpdir / 'yamlmake_logs')
        assert_count_items_dir(tmpdir, 4)
        assert not os.path.exists(tmpdir / 'verybad.txt')

def test_run_until(tmpdir):
    assert_count_items_dir(tmpdir, 0)
    cp = run_yamlmake('basictest.yml', tmpdir, ['--run-until', 'jeremy1'])
    # we instructed it to stop before the task that fails, therefore should return success
    assert cp.returncode == 0
    assert os.path.isfile(tmpdir / 'good0.txt')
    assert os.path.isfile(tmpdir / 'good1.txt')
    assert os.path.isdir(tmpdir / 'yamlmake_logs')
    assert_count_items_dir(tmpdir, 3)

def test_run_from(tmpdir):
    assert_count_items_dir(tmpdir, 0)
    cp = run_yamlmake('basictest.yml', tmpdir, ['--run-from', 'jeremy1'])
    assert cp.returncode == 2
    assert os.path.isfile(tmpdir / 'good1.txt')
    assert os.path.isfile(tmpdir / 'bad.txt')
    assert os.path.isdir(tmpdir / 'yamlmake_logs')
    assert_count_items_dir(tmpdir, 3)


def test_run_only(tmpdir):
    assert_count_items_dir(tmpdir, 0)
    cp = run_yamlmake('basictest.yml', tmpdir, ['--run-only', 'jeremy1'])
    assert cp.returncode == 0
    assert os.path.isfile(tmpdir / 'good1.txt')
    assert os.path.isdir(tmpdir / 'yamlmake_logs')
    assert_count_items_dir(tmpdir, 2)
