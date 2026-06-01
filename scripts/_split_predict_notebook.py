"""One-off: split prediction cells from decomp notebook into a new predict notebook.

Source : notebooks/3. roughness_scale-decomp.ipynb
Target : notebooks/5. roughness_scale-predict.ipynb

Cells moved (by index in source notebook):
    30 (md)  ### OLS regression of n_e/<n> on second-order variance components
    31 (code) OLS table + full vs reduced models
    32 (md)  ### Feature selection: predict n_e/<n> from veg pattern descriptors
    33 (code) Pearson r, mutual info, RF permutation importance
    34 (code) Power-law fit of <S_f>
    37 (md)  ## Tracer flow visuals
    38 (code) format_plot helper
    39 (code) scatter: effect_ratio vs curve, color=aniso
    40 (code) scatter: effect_ratio vs curve, color=sigma
    41 (code) scatter: effect_ratio vs LV, color=fV
    42 (code) dislay_fit helper + predicted-vs-observed plots
    43 (code) groupby l, corr[curve, r_equiv]
    44 (code) groupby l, corr[curve, effect_ratio]

These same indices are deleted from the decomp notebook. Everything else is
untouched.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path("/Users/octaviacrompton/Projects/roughness-scale")
SRC = ROOT / "notebooks" / "3. roughness_scale-decomp.ipynb"
DST = ROOT / "notebooks" / "5. roughness_scale-predict.ipynb"

MOVE_INDICES = [30, 31, 32, 33, 34, 37, 38, 39, 40, 41, 42, 43, 44]
# Setup cells from the decomp notebook to copy verbatim into the new notebook
SETUP_INDICES = [1, 2, 3, 4, 5, 6]


def md_cell(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.splitlines(keepends=True),
    }


def cleared_code_cell(cell: dict) -> dict:
    """Return a copy of a code cell with outputs cleared and exec count reset."""
    c = copy.deepcopy(cell)
    if c.get("cell_type") == "code":
        c["outputs"] = []
        c["execution_count"] = None
    return c


def main() -> None:
    nb = json.loads(SRC.read_text())
    cells = nb["cells"]
    n_before = len(cells)

    # ── 1. Build the new predict notebook ────────────────────────────────────
    title_md = md_cell(
        "## Roughness scaling: predicting $n_e/\\langle n\\rangle$\n"
        "\n"
        "_OLS regression on second-order variance components, feature selection "
        "(Pearson, mutual information, RF permutation importance), power-law fit "
        "of $\\langle S_f\\rangle$, tracer-flow visual diagnostics, and "
        "predicted-vs-observed fits. Extracted from `roughness_scale-decomp.ipynb`._\n"
    )

    # Setup cells: copy verbatim, but rewrite the figure-registry notebook_name
    setup_cells = []
    for i in SETUP_INDICES:
        c = cleared_code_cell(cells[i])
        if c.get("cell_type") == "code":
            src = "".join(c["source"])
            src = src.replace(
                "notebook_name='roughness_scale-decomp.ipynb'",
                "notebook_name='5. roughness_scale-predict.ipynb'",
            )
            c["source"] = src.splitlines(keepends=True)
        setup_cells.append(c)

    moved_cells = [cleared_code_cell(cells[i]) for i in MOVE_INDICES]

    new_nb = {
        "cells": [title_md, *setup_cells, *moved_cells],
        "metadata": copy.deepcopy(nb.get("metadata", {})),
        "nbformat": nb.get("nbformat", 4),
        "nbformat_minor": nb.get("nbformat_minor", 5),
    }

    DST.write_text(json.dumps(new_nb, indent=1, ensure_ascii=False) + "\n")
    print(f"Wrote {DST}  ({len(new_nb['cells'])} cells)")

    # ── 2. Delete moved cells from the decomp notebook ───────────────────────
    keep = [c for i, c in enumerate(cells) if i not in set(MOVE_INDICES)]
    nb["cells"] = keep
    SRC.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
    print(f"Updated {SRC}: {n_before} → {len(keep)} cells (removed {len(MOVE_INDICES)})")


if __name__ == "__main__":
    main()
