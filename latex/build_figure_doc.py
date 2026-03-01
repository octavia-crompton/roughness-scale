#!/usr/bin/env python3
"""
Generate a LaTeX figure document from a figure_registry_concise.txt file.

Usage:
    python build_figure_doc.py              # builds all batch folders
    python build_figure_doc.py runaround_smooth   # build for one batch

The script:
  1. Parses figures/<batch>/figure_registry_concise.txt
  2. Writes   latex/<batch>_figures.tex
  3. Compiles with pdflatex (if available)

Called automatically by watch_registry.sh when the concise registry changes.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # roughness-scale/
FIGURES_DIR = ROOT / "figures"
LATEX_DIR   = ROOT / "latex"

# MacTeX on macOS installs here; add to PATH so pdflatex is found.
_TEXBIN = "/Library/TeX/texbin"
if os.path.isdir(_TEXBIN) and _TEXBIN not in os.environ.get("PATH", ""):
    os.environ["PATH"] = _TEXBIN + ":" + os.environ.get("PATH", "")

# ── Registry parser ────────────────────────────────────────────────────────────
_HEADER_RE = re.compile(
    r"^(Fig\s+\d+|SI\s*\d+)\s*—\s*(.+?)$",
    re.MULTILINE,
)


def parse_concise_registry(registry_path: Path) -> list[dict]:
    """Return list of {fig_id, filename, caption} dicts, in file order."""
    text = registry_path.read_text()
    entries = []
    headers = list(_HEADER_RE.finditer(text))
    for i, m in enumerate(headers):
        fig_id  = m.group(1).strip()          # "Fig 4", "SI 1"
        fname   = m.group(2).strip()          # "fig4_obs_vs_pred_re_6panel.png"
        # Caption is the indented text between this header and the next (or EOF)
        start = m.end()
        end   = headers[i + 1].start() if i + 1 < len(headers) else len(text)
        block = text[start:end].strip()
        caption = block if block else fname   # fallback to filename if caption empty
        entries.append(dict(fig_id=fig_id, filename=fname, caption=caption))
    return entries


# ── LaTeX generation ───────────────────────────────────────────────────────────
def _latex_escape(s: str) -> str:
    """Escape special LaTeX characters in plain text, preserving $..$ math."""
    # Split on $…$ to keep inline math intact
    parts = re.split(r'(\$[^$]+\$)', s)
    escaped = []
    for i, part in enumerate(parts):
        if i % 2 == 1:  # inside $…$
            escaped.append(part)
        else:
            part = part.replace('\\', r'\textbackslash{}')
            for ch in ['&', '%', '#', '_', '{', '}']:
                part = part.replace(ch, f'\\{ch}')
            part = part.replace('~', r'\textasciitilde{}')
            part = part.replace('^', r'\textasciicircum{}')
            # Handle ≈ → \approx  (common in our captions)
            part = part.replace('≈', r'$\approx$')
            escaped.append(part)
    return ''.join(escaped)


def _fig_label(fig_id: str) -> str:
    """'Fig 4' → 'fig:fig4', 'SI 1' → 'fig:SI1'."""
    return "fig:" + fig_id.replace(" ", "").lower().replace("fig", "fig")


def generate_tex(batch: str, entries: list[dict]) -> str:
    """Build the full .tex source for one batch."""
    rel_fig_dir = f"../figures/{batch}"
    body_lines = []

    for e in entries:
        label   = _fig_label(e["fig_id"])
        fname   = e["filename"]
        caption = _latex_escape(e["caption"])
        body_lines.append(rf"""
\begin{{figure}}[htbp]
  \centering
  \includegraphics[width=\textwidth]{{{rel_fig_dir}/{fname}}}
  \caption{{{caption}}}
  \label{{{label}}}
\end{{figure}}
""")

    figures_block = "\n".join(body_lines)
    title_escaped = _latex_escape(batch.replace("_", " ").title())

    tex = rf"""\documentclass[11pt]{{article}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{hyperref}}
\usepackage[font=small,labelfont=bf]{{caption}}

\title{{Figure Compilation — {title_escaped}}}
\author{{Auto-generated from figure\_registry\_concise.txt}}
\date{{\today}}

\begin{{document}}
\maketitle
\listoffigures
\clearpage
{figures_block}
\end{{document}}
"""
    return tex


# ── Compile ────────────────────────────────────────────────────────────────────
def compile_tex(tex_path: Path) -> bool:
    """Run pdflatex twice (for TOC). Returns True on success.
    Silently skips if pdflatex is not installed.
    """
    import shutil
    if shutil.which("pdflatex") is None:
        print("  ⚠ pdflatex not found — .tex written but not compiled.")
        print("    Install a TeX distribution:  brew install --cask mactex-no-gui")
        return False
    cmd = [
        "pdflatex",
        "-interaction=nonstopmode",
        "-halt-on-error",
        "-output-directory", str(tex_path.parent),
        str(tex_path),
    ]
    for pass_num in (1, 2):
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tex_path.parent))
        if result.returncode != 0:
            print(f"  ✗ pdflatex pass {pass_num} failed for {tex_path.name}")
            print(result.stdout[-1500:] if len(result.stdout) > 1500 else result.stdout)
            return False
    print(f"  ✓ {tex_path.with_suffix('.pdf').name}")
    return True


# ── Main ───────────────────────────────────────────────────────────────────────
def build_batch(batch: str) -> bool:
    """Parse registry, write .tex, compile .pdf for one batch folder."""
    registry = FIGURES_DIR / batch / "figure_registry_concise.txt"
    if not registry.exists():
        print(f"  ⚠ No concise registry in figures/{batch}/ — skipping")
        return False

    entries = parse_concise_registry(registry)
    if not entries:
        print(f"  ⚠ No figure entries parsed from {registry} — skipping")
        return False

    tex_src  = generate_tex(batch, entries)
    tex_path = LATEX_DIR / f"{batch}_figures.tex"
    tex_path.write_text(tex_src)
    print(f"  Wrote {tex_path.relative_to(ROOT)}")

    return compile_tex(tex_path)


def main():
    LATEX_DIR.mkdir(exist_ok=True)

    if len(sys.argv) > 1:
        batches = sys.argv[1:]
    else:
        # auto-discover all batch folders that contain a concise registry
        batches = sorted(
            d.name for d in FIGURES_DIR.iterdir()
            if d.is_dir() and (d / "figure_registry_concise.txt").exists()
        )

    if not batches:
        print("No figure batches found.")
        return

    for batch in batches:
        print(f"\n── {batch} ──")
        build_batch(batch)


if __name__ == "__main__":
    main()
