"""Manifests store absolute paths; loaders must survive a moved or renamed result folder."""
from pathlib import Path

from fastccc.results import _relocate, _explorer_path_from_manifest


def test_relocate_prefers_existing_stored_path(tmp_path):
    stored = tmp_path / 'stored'
    stored.mkdir()
    assert _relocate(str(stored), tmp_path / 'other') == stored


def test_relocate_falls_back_to_candidates(tmp_path):
    gone = tmp_path / 'old_name' / 'groupA'
    moved = tmp_path / 'new_name' / 'groupA'
    moved.mkdir(parents=True)
    assert _relocate(str(gone), tmp_path / 'new_name' / gone.name) == moved


def test_relocate_returns_stored_when_nothing_exists(tmp_path):
    gone = tmp_path / 'old' / 'x'
    assert _relocate(str(gone), tmp_path / 'also_missing') == gone
    assert _relocate(None, tmp_path / 'also_missing') is None


def test_explorer_path_relocates_into_result_dir(tmp_path):
    result_dir = tmp_path / 'renamed'
    result_dir.mkdir()
    (result_dir / 'explorer.html').write_text('<html></html>')
    manifest = {'outputs': {'explorer_html': str(tmp_path / 'original' / 'explorer.html')}}
    assert _explorer_path_from_manifest(manifest, result_dir) == result_dir / 'explorer.html'
    assert _explorer_path_from_manifest({'outputs': {}}, result_dir) is None


def test_explore_records_path_without_manifest(tmp_path, monkeypatch):
    """``res.explore()`` must work on a result folder without a manifest (the property has no setter)."""
    from fastccc import explorer as explorer_mod
    from fastccc.results import RunResult
    out = tmp_path / 'explorer.html'
    monkeypatch.setattr(explorer_mod, 'build_explorer', lambda result_dir, output=None, **kw: out)
    res = RunResult(result_dir=tmp_path, task_id='abc', name='x')
    assert res.explorer_path is None
    assert res.explore() == out
    assert res.explorer_path == out


def test_legacy_database_paths_resolve_to_packaged_copy(tmp_path):
    """Original entry points accept the former checkout paths and plain names."""
    from fastccc.resources import legacy_database_path, PACKAGED_DB_DIR
    for raw in ('./db/CPDBv4.1.0/', './db/CPDBv5.0.0', 'CPDBv5.0.0'):
        name = raw.rstrip('/').split('/')[-1]
        assert legacy_database_path(raw) == str((PACKAGED_DB_DIR / name).resolve())
    own = tmp_path / 'mydb'; own.mkdir(); (own / 'interaction_table.csv').write_text('x\n')
    assert legacy_database_path(str(own) + '/') == str(own.resolve())          # an existing folder is used as given
    import inspect, fastccc
    for f in (fastccc.statistical_analysis_method, fastccc.Cauchy_combination_of_statistical_analysis_methods,
              fastccc.build_reference.build_reference_workflow, fastccc.infer_query.infer_query_workflow):
        assert 'legacy_database_path(database_file_path)' in inspect.getsource(f)
