"""
Figure registry helpers for the roughness-scale notebooks.

Stores figure metadata (filename, description, timestamp) in two plain-text
files alongside the saved PNGs:
  - figure_registry.txt         full, detailed entries
  - figure_registry_concise.txt short two-sentence summaries

Usage (in a notebook)
---------------------
    import sys
    sys.path.insert(0, "/Users/octaviacrompton/Projects/roughness-scale/src")
    import figure_registry as fig_reg

    # call once after out_dir is set
    fig_reg.configure(out_dir, notebook_name='roughness_scale-analysis.ipynb')

    # re-expose convenience names so existing cells don't change
    _fig_dirs             = fig_reg._fig_dirs
    update_figure_registry = fig_reg.update_figure_registry
"""

import os as _os
import re as _re
from datetime import datetime as _datetime

# ── module-level state (set via configure()) ──────────────────────────────────
_OUT_DIR       : str = ""
_NOTEBOOK_NAME : str = "roughness_scale-analysis.ipynb"


def configure(out_dir: str, notebook_name: str = _NOTEBOOK_NAME) -> None:
    """Initialise the registry for a session.

    Call once after ``out_dir`` is known (typically right after loading data).

    Parameters
    ----------
    out_dir : str
        Path to the simulation batch directory (e.g. ``.../Tests/runaround_smooth``).
    notebook_name : str, optional
        Filename of the source notebook (default: ``roughness_scale-analysis.ipynb``).
    """
    global _OUT_DIR, _NOTEBOOK_NAME
    _OUT_DIR       = out_dir
    _NOTEBOOK_NAME = notebook_name


# ── directory helpers ─────────────────────────────────────────────────────────

def _fig_dirs():
    """Return ``(fig_dir, scratch_dir, registry_path)`` derived from *_OUT_DIR*."""
    if not _OUT_DIR:
        raise RuntimeError(
            "figure_registry.configure(out_dir) must be called before _fig_dirs()."
        )
    batch       = _os.path.basename(_OUT_DIR.rstrip("/"))
    fig_dir     = _os.path.join("..", "figures", batch)
    scratch_dir = _os.path.join("..", "figures", batch, "scratch")
    registry    = _os.path.join(fig_dir, "figure_registry.txt")
    _os.makedirs(fig_dir,     exist_ok=True)
    _os.makedirs(scratch_dir, exist_ok=True)
    return fig_dir, scratch_dir, registry


# ── sort & parse helpers ──────────────────────────────────────────────────────

def _sort_key(fig_id: str):
    """Sort key: main figures first (fig1 … figN), then SI figures."""
    lo   = fig_id.lower()
    is_si = lo.startswith("si") or lo.startswith("figsi") or lo.startswith("fig_si")
    m    = _re.search(r'(\d+)', fig_id)
    num  = int(m.group(1)) if m else 999
    return (1 if is_si else 0, num, fig_id)


def _parse_registry(registry_path: str) -> dict:
    """Read existing registry entries.  Returns ``{fig_id: entry_dict}``."""
    entries = {}
    if not _os.path.exists(registry_path):
        return entries
    with open(registry_path) as f:
        text = f.read()

    # ── current 5-field header ────────────────────────────────────────────
    for m in _re.finditer(
        r"^### (\S+) ###\n"
        r"File     : (.+)\n"
        r"Updated  : (.+)\n"
        r"Notebook : (.+)\n"
        r"Saved in : (.+)\n"
        r"─+\n"
        r"(.*?)"
        r"^### end \1 ###",
        text, _re.DOTALL | _re.MULTILINE,
    ):
        fid = m.group(1)
        entries[fid] = dict(
            filename  = m.group(2).strip(),
            updated   = m.group(3).strip(),
            notebook  = m.group(4).strip().replace("notebooks/", ""),
            save_dir  = m.group(5).strip(),
            description = m.group(6).strip(),
            concise   = "",
        )

    # ── legacy 2-field header fallback ───────────────────────────────────
    if not entries:
        for m in _re.finditer(
            r"^### (\S+) ###\n"
            r"File    : (.+)\n"
            r"Updated : (.+)\n"
            r"─+\n"
            r"(.*?)"
            r"^### end \1 ###",
            text, _re.DOTALL | _re.MULTILINE,
        ):
            fid = m.group(1)
            entries[fid] = dict(
                filename    = m.group(2).strip(),
                updated     = m.group(3).strip(),
                description = m.group(4).strip(),
                concise     = "",
            )

    # ── back-fill concise text ────────────────────────────────────────────
    concise_path = registry_path.replace("figure_registry.txt",
                                          "figure_registry_concise.txt")
    if _os.path.exists(concise_path):
        with open(concise_path) as f:
            ctext = f.read()
        _CONC_RE = _re.compile(
            r"^(Fig\s+\d+|SI\s*\d+)\s*—\s*.+?\n\s{2}(.+?)$",
            _re.MULTILINE,
        )
        for cm in _CONC_RE.finditer(ctext):
            tag     = cm.group(1).strip()
            caption = cm.group(2).strip()
            fid     = tag.replace(" ", "").replace("Fig", "fig")
            if fid in entries and not entries[fid]["concise"]:
                entries[fid]["concise"] = caption

    return entries


# ── write helpers ─────────────────────────────────────────────────────────────

def _rewrite_registry(entries: dict, registry_path: str) -> None:
    """Write both ``figure_registry.txt`` and ``figure_registry_concise.txt``."""
    divider    = "─" * 72
    sorted_ids = sorted(entries, key=_sort_key)
    fig_dir, _, _ = _fig_dirs()

    # ── full registry ─────────────────────────────────────────────────────
    lines = []
    lines.append(
        f"Figure registry  •  created {_datetime.now().strftime('%Y-%m-%d %H:%M')}"
    )
    lines.append(f"Source notebook  :  notebooks/{_NOTEBOOK_NAME}")
    lines.append(f"Figures saved in :  {_os.path.abspath(fig_dir)}")
    lines.append(divider)
    lines.append("")
    for fid in sorted_ids:
        e = entries[fid]
        lines.append(f"### {fid} ###")
        lines.append(f"File     : {e['filename']}")
        lines.append(f"Updated  : {e['updated']}")
        lines.append(f"Notebook : notebooks/{e.get('notebook', _NOTEBOOK_NAME)}")
        lines.append(f"Saved in : {e.get('save_dir', fig_dir)}")
        lines.append(divider)
        lines.append(e["description"].strip())
        lines.append(f"### end {fid} ###")
        lines.append("")
    with open(registry_path, "w") as f:
        f.write("\n".join(lines))

    # ── concise registry ──────────────────────────────────────────────────
    concise_path = registry_path.replace("figure_registry.txt",
                                          "figure_registry_concise.txt")
    clines = []
    clines.append(
        f"Figure Registry (concise)  •  {_datetime.now().strftime('%Y-%m-%d %H:%M')}"
    )
    clines.append(
        f"Source: notebooks/{_NOTEBOOK_NAME}  |  "
        f"Figures: {_os.path.abspath(fig_dir)}"
    )
    clines.append("=" * 72)
    clines.append("")
    for fid in sorted_ids:
        e   = entries[fid]
        tag = fid.upper().replace("FIG", "Fig ")
        short = e.get("concise", "") or e["description"].split("\n")[0]
        clines.append(f"{tag} — {e['filename']}")
        clines.append(f"  {short}")
        clines.append("")
    with open(concise_path, "w") as f:
        f.write("\n".join(clines))


# ── public API ────────────────────────────────────────────────────────────────

def update_figure_registry(
    fig_id: str,
    filename: str,
    description: str,
    concise: str = "",
) -> None:
    """Write or replace the section for *fig_id* in both registries.

    Parameters
    ----------
    fig_id : str
        Identifier such as ``'fig1'``, ``'fig4'``, ``'SI1'``.
    filename : str
        Name of the saved PNG/PDF file (basename only).
    description : str
        Full multi-line description for the main registry.
    concise : str, optional
        Two-sentence human summary for the concise registry.
        Falls back to the first line of *description*.
    """
    _, _, registry_path = _fig_dirs()
    now     = _datetime.now().strftime("%Y-%m-%d %H:%M")
    entries = _parse_registry(registry_path)
    entries[fig_id] = dict(
        filename    = filename,
        updated     = now,
        description = description.strip(),
        concise     = concise.strip() if concise else description.strip().split("\n")[0],
    )
    _rewrite_registry(entries, registry_path)
    print(f"Registry updated → {registry_path}  [{fig_id}]")
