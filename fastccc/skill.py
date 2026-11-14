"""Locate, show and install the FastCCC agent skill that ships inside the package.

The skill (``SKILL.md`` plus ``references/``) teaches coding agents such as Claude
Code, Codex or Cursor how to run FastCCC from plain-language requests.  It is
also available with ``npx skills add Svvord/FastCCC``.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Iterable, Optional

SKILL_NAME = 'fastccc'
PACKAGED_SKILL_DIR = Path(__file__).resolve().parent / 'data' / 'skill' / SKILL_NAME

# agent name -> (project-level skills dir, user-level skills dir)
AGENT_SKILL_DIRS: dict[str, tuple[str, str]] = {
    'claude': ('.claude/skills', '~/.claude/skills'),
    'codex': ('.agents/skills', '~/.agents/skills'),
    'cursor': ('.agents/skills', '~/.cursor/skills'),
    'opencode': ('.agents/skills', '~/.config/opencode/skills'),
    'gemini': ('.agents/skills', '~/.gemini/skills'),
    'copilot': ('.agents/skills', '~/.copilot/skills'),
}
DEFAULT_AGENTS = ('claude', 'codex')


def skill_path() -> Path:
    path = PACKAGED_SKILL_DIR / 'SKILL.md'
    if not path.is_file():
        raise FileNotFoundError(f"packaged skill not found at {path}")
    return path


def skill_text() -> str:
    return skill_path().read_text(encoding='utf-8')


def skill_targets(agents: Iterable[str] = DEFAULT_AGENTS, global_: bool = False,
                  project_dir: Optional[str | os.PathLike] = None) -> dict[Path, list[str]]:
    """Resolve the skill directories for the requested agents (deduplicated)."""
    names = list(agents)
    if names == ['all'] or 'all' in names:
        names = list(AGENT_SKILL_DIRS)
    unknown = [a for a in names if a not in AGENT_SKILL_DIRS]
    if unknown:
        raise ValueError(f"unknown agent(s) {unknown}; choose from {', '.join(AGENT_SKILL_DIRS)} or 'all'")
    base = Path(project_dir).expanduser().resolve() if project_dir else Path.cwd()
    targets: dict[Path, list[str]] = {}
    for a in names:
        project_rel, global_abs = AGENT_SKILL_DIRS[a]
        root = Path(global_abs).expanduser() if global_ else base / project_rel
        targets.setdefault(root / SKILL_NAME, []).append(a)
    return targets


def install_skill(agents: Iterable[str] = DEFAULT_AGENTS, global_: bool = False,
                  project_dir: Optional[str | os.PathLike] = None, force: bool = False) -> list[dict]:
    """Copy the packaged skill into the agents' skill folders.

    Returns one record per target: ``{'path', 'agents', 'status'}`` with status
    ``installed``, ``updated`` or ``exists`` (left untouched; pass ``force`` to overwrite).
    """
    src = skill_path().parent
    records = []
    for target, names in skill_targets(agents, global_, project_dir).items():
        existed = (target / 'SKILL.md').is_file()
        if existed and not force:
            records.append({'path': str(target), 'agents': names, 'status': 'exists'})
            continue
        if existed:
            shutil.rmtree(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src, target)
        records.append({'path': str(target), 'agents': names, 'status': 'updated' if existed else 'installed'})
    return records
