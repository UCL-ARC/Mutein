import yamlmake as ym
import pytest
import os
import datetime
from datetime import timezone
from collections import namedtuple

FakeArgs = namedtuple('FakeArgs', 'yaml log_dir run_only run_from run_until dry_run quiet no_logs')

def test_overall():
    assert 3 == 3

@pytest.mark.parametrize("num_reps", [1, 2])
def test_basic_workflow(num_reps):
    yamlfile = os.path.join(os.path.dirname(__file__), 'basictest.yml')
    pipeline = ym.parse_yaml(yamlfile)
    for it in range(num_reps):
        # run the workflow more than once to check that if output files were created by a
        # failed command, on subsequent attempts to run the workflow those files are not used and the
        # workflow continues to fail
        args = FakeArgs(None, None, None, None, None, None, None, None)
        ym.process(pipeline, yamlfile, args=args)
        assert os.path.isfile('good.txt')
        assert os.path.isfile('bad.txt')
        o = os.path.getmtime('bad.txt')
        e = datetime.datetime(1980, 1, 1, tzinfo=timezone.utc).timestamp()
        assert e == o
        assert not os.path.exists('verybad.txt')