"""`fastccc run` module-selection flags."""
from typer.testing import CliRunner

from fastccc.cli import app


def test_run_help_lists_module_flags():
    out = CliRunner().invoke(app, ['run', '--help']).output
    for flag in ('--summary', '--aggregation', '--combination'):
        assert flag in out


def test_run_rejects_unknown_summary(tmp_path):
    result = CliRunner().invoke(app, ['run', str(tmp_path / 'missing.h5ad'), '-c', 'cell_type', '--summary', 'Mean', '--summary', 'Mode'])
    assert result.exit_code != 0
    assert 'Mode' in result.output and 'Quantile_<q>' in result.output


def test_run_accepts_custom_quantile_summary(tmp_path):
    result = CliRunner().invoke(app, ['run', str(tmp_path / 'missing.h5ad'), '-c', 'cell_type', '--summary', 'Quantile_0.8', '--summary', 'Q2'])
    assert 'not in' not in result.output  # validation passed; the run then fails on the missing file
    assert result.exit_code != 0


def test_cli_json_output_is_clean_stdout():
    """Progress logs go to stderr so that ``--json`` on stdout parses as a whole."""
    import json
    import subprocess
    import sys
    r = subprocess.run([sys.executable, '-m', 'fastccc.cli', 'db', 'list', '--json'], capture_output=True, text=True, check=True)
    data = json.loads(r.stdout)
    assert any(d['name'] == 'CPDBv5.0.0' for d in data['databases'])
