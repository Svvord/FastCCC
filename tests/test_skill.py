"""The packaged agent skill must stay in sync with the CLI/API and with the repo copies."""
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

import fastccc
from fastccc import skill as sk
from fastccc.cli import app

REPO = Path(__file__).resolve().parents[1]
KNOWN_COMMANDS = {'inspect', 'run', 'ref', 'build-ref', 'map', 'panels', 'report', 'explore', 'skill', 'db'}


def _frontmatter(text: str) -> dict:
    m = re.match(r'^---\n(.*?)\n---\n', text, flags=re.S)
    assert m, 'SKILL.md must start with a YAML front matter block'
    out = {}
    for line in m.group(1).splitlines():
        key, _, value = line.partition(':')
        out[key.strip()] = value.strip()
    return out


def test_frontmatter_and_name():
    fm = _frontmatter(sk.skill_text())
    assert fm['name'] == sk.SKILL_NAME == 'fastccc'
    assert 100 < len(fm['description']) <= 1024
    for word in ('cell-cell communication', 'ligand-receptor', 'healthy tissue', 'Up or Down'):
        assert word in fm['description']


def _skill_files(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): p.read_text(encoding='utf-8')
            for p in sorted(root.rglob('*')) if p.is_file() and p.name != 'openai.yaml'}


@pytest.mark.skipif(not (REPO / '.claude' / 'skills').is_dir(), reason='repo checkout only')
def test_repo_copies_match_packaged_skill():
    packaged = _skill_files(sk.PACKAGED_SKILL_DIR)
    assert 'SKILL.md' in packaged and 'references/workflows.md' in packaged
    for copy in (REPO / '.claude' / 'skills' / 'fastccc', REPO / '.agents' / 'skills' / 'fastccc'):
        assert _skill_files(copy) == packaged, f'{copy} differs from the packaged skill'
    assert (REPO / '.agents' / 'skills' / 'fastccc' / 'agents' / 'openai.yaml').is_file()


def test_every_mentioned_cli_command_exists():
    text = sk.skill_text() + (sk.PACKAGED_SKILL_DIR / 'references' / 'workflows.md').read_text(encoding='utf-8')
    mentioned = set(re.findall(r'\bfastccc ([a-z][a-z-]*)\b', text))
    unknown = mentioned - KNOWN_COMMANDS
    assert not unknown, f'skill mentions unknown fastccc subcommands: {unknown}'
    runner = CliRunner()
    for cmd in mentioned:
        assert runner.invoke(app, [cmd, '--help']).exit_code == 0, cmd
    options = set(re.findall(r'(--[a-z][a-z-]+)', text))
    helps = {cmd: runner.invoke(app, [cmd, '--help']).output for cmd in KNOWN_COMMANDS}
    helps['skill install'] = runner.invoke(app, ['skill', 'install', '--help']).output
    for sub in ('list', 'show', 'validate', 'import'):
        helps[f'db {sub}'] = runner.invoke(app, ['db', sub, '--help']).output
    helps['top'] = runner.invoke(app, ['--help']).output
    for opt in options:
        assert any(opt in h for h in helps.values()), f'option {opt} mentioned in the skill is not a CLI option'


def test_every_mentioned_python_name_exists():
    text = sk.skill_text() + (sk.PACKAGED_SKILL_DIR / 'references' / 'workflows.md').read_text(encoding='utf-8')
    names = set(re.findall(r'\bfastccc\.([A-Za-z_]+)', text))
    for name in names:
        assert hasattr(fastccc, name), f'fastccc.{name} is mentioned in the skill but does not exist'
    def has(cls, attr):
        return hasattr(cls, attr) or attr in getattr(cls, '__dataclass_fields__', {})
    for attr in ('trend_counts', 'interactions', 'interaction_summary', 'celltype_mapping', 'summary', 'report', 'reference', 'results', 'report_path'):
        assert has(fastccc.ReferenceResult, attr), attr
    for attr in ('summary', 'significant', 'pvals', 'strength', 'report', 'report_path'):
        assert has(fastccc.RunResult, attr), attr
    for attr in ('groups', 'all', 'report_path'):
        assert has(fastccc.MultiRunResult, attr), attr
    for attr in ('format', 'ambiguous', 'unmatched_query', 'mapping'):
        assert has(fastccc.MappingSuggestion, attr), attr
    for attr in ('celltypes', 'n_cells', 'path'):
        assert has(fastccc.ReferencePanel, attr), attr


def test_install_skill_into_project(tmp_path):
    recs = sk.install_skill(project_dir=tmp_path)
    paths = {Path(r['path']) for r in recs}
    assert paths == {tmp_path / '.claude' / 'skills' / 'fastccc', tmp_path / '.agents' / 'skills' / 'fastccc'}
    assert all(r['status'] == 'installed' for r in recs)
    assert (tmp_path / '.claude' / 'skills' / 'fastccc' / 'references' / 'workflows.md').is_file()
    again = sk.install_skill(project_dir=tmp_path)
    assert all(r['status'] == 'exists' for r in again)
    forced = sk.install_skill(agents=['cursor', 'codex'], project_dir=tmp_path, force=True)
    assert len(forced) == 1 and forced[0]['status'] == 'updated' and set(forced[0]['agents']) == {'cursor', 'codex'}
    with pytest.raises(ValueError, match='unknown agent'):
        sk.install_skill(agents=['vim'], project_dir=tmp_path)
    assert len(sk.skill_targets(['all'], project_dir=tmp_path)) == 2


def test_cli_skill_commands(tmp_path):
    runner = CliRunner()
    out = runner.invoke(app, ['skill', 'show'])
    assert out.exit_code == 0 and out.output.startswith('---\nname: fastccc')
    out = runner.invoke(app, ['skill', 'install', '--dir', str(tmp_path), '-a', 'claude'])
    assert out.exit_code == 0 and (tmp_path / '.claude' / 'skills' / 'fastccc' / 'SKILL.md').is_file(), out.output
    assert 'npx skills add Svvord/FastCCC' in out.output
