"""Database tooling: build from pairs, validate, metadata, import and name resolution."""
import json

import pandas as pd
import pytest
from typer.testing import CliRunner

from fastccc import dbtools, resources
from fastccc.cli import app

PAIRS = pd.DataFrame({
    'ligand': ['TGFB1', 'TGFB1', 'IL6', 'CLCF1_CRLF1', 'Dll1'],
    'receptor': ['TGFBR1_TGFBR2', 'TGFBR1_TGFBR2', 'IL6R_IL6ST', 'CNTFR_IL6ST', 'Notch1'],
    'classification': ['TGF-beta', 'TGF-beta', 'IL6 family', None, None],
})


def test_split_subunits():
    assert dbtools.split_subunits('ITGA1_ITGB1') == ('ITGA1', 'ITGB1')
    assert dbtools.split_subunits('COMPLEX:P05556_P53708') == ('P05556', 'P53708')
    assert dbtools.split_subunits('A, B') == ('A', 'B')
    assert dbtools.split_subunits(' ') == ()


def test_build_validate_and_info(tmp_path):
    out = dbtools.build_database(PAIRS, 'toy', tmp_path, species='human', source='unit test', version='0')
    rep = dbtools.validate_database(out)
    assert rep['ok'], rep
    s = rep['summary']
    assert s['n_interactions'] == 4           # the duplicated TGFB1 pair collapses
    assert s['n_complexes'] == 4              # TGFBR1_TGFBR2, IL6R_IL6ST, CLCF1_CRLF1, CNTFR_IL6ST
    assert s['has_classification'] is True
    inter = pd.read_csv(out / 'interaction_table.csv')
    genes = pd.read_csv(out / 'gene_table.csv')
    assert genes.loc[genes.hgnc_symbol == 'TGFB1', 'ensembl'].item().startswith('ENSG')   # human mapping filled
    assert set(inter.columns) >= {'id_cp_interaction', 'multidata_1_id', 'multidata_2_id', 'classification'}
    meta = dbtools.database_info(out)
    assert meta['species'] == 'human' and meta['source'] == 'unit test' and meta['n_interactions'] == 4
    # the core loader understands the layout
    from fastccc.preprocess import get_interactions
    assert len(get_interactions(str(out))) == 4


def test_build_rejects_bad_input(tmp_path):
    with pytest.raises(ValueError):
        dbtools.build_database(pd.DataFrame({'ligand': ['A']}), 'x', tmp_path)
    with pytest.raises(ValueError):
        dbtools.build_database(PAIRS, 'x', tmp_path, species='martian')
    dbtools.build_database(PAIRS, 'x', tmp_path)
    with pytest.raises(FileExistsError):
        dbtools.build_database(PAIRS, 'x', tmp_path)


def test_validate_reports_broken_references(tmp_path):
    out = dbtools.build_database(PAIRS, 'toy', tmp_path)
    inter = pd.read_csv(out / 'interaction_table.csv')
    inter.loc[0, 'multidata_2_id'] = 9999
    inter.to_csv(out / 'interaction_table.csv', index=False)
    rep = dbtools.validate_database(out)
    assert not rep['ok'] and any('neither a protein nor a complex' in e for e in rep['errors'])
    assert not dbtools.validate_database(tmp_path / 'missing')['ok']


@pytest.mark.parametrize('name', resources.list_databases(include_user=False))
def test_packaged_databases_validate(name):
    rep = dbtools.validate_database(resources.PACKAGED_DB_DIR / name)
    assert rep['ok'], rep['errors']
    meta = dbtools.database_info(resources.PACKAGED_DB_DIR / name)
    assert meta['species'] in dbtools.SPECIES and meta.get('citation')


def test_import_and_resolve_user_database(tmp_path, monkeypatch):
    monkeypatch.setenv('FASTCCC_HOME', str(tmp_path / 'home'))
    table = tmp_path / 'pairs.tsv'
    PAIRS.rename(columns={'ligand': 'L', 'receptor': 'R', 'classification': 'pathway'}).to_csv(table, sep='\t', index=False)
    path = dbtools.import_database(table, name='mine', species='mouse',
                                   pairs_kwargs=dict(ligand_col='L', receptor_col='R', classification_col='pathway'))
    assert path == (tmp_path / 'home' / 'databases' / 'mine').resolve()
    assert 'mine' in resources.list_databases()
    assert resources.resolve_database('mine') == path
    assert resources.resolve_database('./somewhere/mine') == path     # last path component falls back to the name
    assert dbtools.database_species(path) == 'mouse'
    # a directory import copies and re-validates
    copied = dbtools.import_database(path, name='mine2', citation='x')
    assert dbtools.database_info(copied)['citation'] == 'x'
    with pytest.raises(FileExistsError):
        dbtools.import_database(path, name='mine2')
    with pytest.raises(FileNotFoundError):
        resources.resolve_database('nothing-here')


def test_cli_db_commands(tmp_path, monkeypatch):
    monkeypatch.setenv('FASTCCC_HOME', str(tmp_path / 'home'))
    runner = CliRunner()
    out = runner.invoke(app, ['db', 'list', '--json']).output
    names = [d['name'] for d in json.loads(out)['databases']]
    assert 'CPDBv5.0.0' in names
    out = runner.invoke(app, ['db', 'show', 'CPDBv5.0.0']).output
    assert 'species' in out and 'citation' in out
    table = tmp_path / 'pairs.csv'
    PAIRS.to_csv(table, index=False)
    r = runner.invoke(app, ['db', 'import', str(table), '--name', 'toy', '--classification-col', 'classification', '--source', 'test'])
    assert r.exit_code == 0, r.output
    assert 'toy' in runner.invoke(app, ['db', 'list']).output
    r = runner.invoke(app, ['db', 'validate', str(tmp_path / 'home' / 'databases' / 'toy')])
    assert r.exit_code == 0 and 'OK' in r.output
    r = runner.invoke(app, ['db', 'validate', str(tmp_path)])
    assert r.exit_code == 1
    r = runner.invoke(app, ['db', 'import', str(table), '--name', 'toy'])
    assert r.exit_code == 1 and 'exists' in r.output


def test_database_fit_check_rejects_wrong_species():
    import anndata
    import numpy as np
    import pytest
    from fastccc.api import _check_database_fit
    from fastccc.resources import resolve_database
    mouse_like = anndata.AnnData(np.ones((5, 4), dtype=np.float32),
                                 var=pd.DataFrame(index=['Il2', 'Il2ra', 'Tgfb1', 'Gm1234']))
    with pytest.raises(ValueError, match='gene'):
        _check_database_fit(resolve_database('CPDBv5.0.0'), mouse_like)
    db = resolve_database('CellChatv1.6.1')
    some = pd.read_csv(db / 'gene_table.csv')['hgnc_symbol'].drop_duplicates().head(40).tolist()
    human = anndata.AnnData(np.ones((5, len(some)), dtype=np.float32), var=pd.DataFrame(index=some))
    fit = _check_database_fit(db, human)
    assert 0.05 <= fit['gene_coverage'] < 0.5 and fit['species'] == 'human'   # small panel: warns, does not raise


def test_database_fit_accepts_ensembl_ids():
    import anndata
    import numpy as np
    from fastccc.api import _check_database_fit
    from fastccc.resources import resolve_database
    db = resolve_database('CPDBv5.0.0')
    ens = pd.read_csv(db / 'gene_table.csv')['ensembl'].dropna().astype(str).drop_duplicates().head(200).tolist()
    data = anndata.AnnData(np.ones((5, len(ens)), dtype=np.float32), var=pd.DataFrame(index=ens))
    with pytest.raises(ValueError, match="convert_type='ensembl'"):
        _check_database_fit(db, data)                     # symbol mode: no overlap, but the hint names the fix
    fit = _check_database_fit(db, data, 'ensembl')
    assert fit['gene_id'] == 'ensembl' and fit['n_database_genes_found'] == 200
    assert 0.1 < dbtools.gene_coverage(db, ens) < 0.2    # auto-detects Ensembl ids


def test_cli_run_gene_id_option(tmp_path):
    runner = CliRunner()
    assert '--gene-id' in runner.invoke(app, ['run', '--help']).output
    r = runner.invoke(app, ['run', str(tmp_path / 'missing.h5ad'), '-c', 'cell_type', '--gene-id', 'entrez'])
    assert r.exit_code == 1 and 'ensembl' in r.output
