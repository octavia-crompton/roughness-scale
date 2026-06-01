#!/usr/bin/env python3
"""Patch 4. roughness_scale-cf-channel_eqs.ipynb:
1. Fix Yen titles: wrap R^{1/6} and R^{1/3} in LaTeX $...$
2. Increase display DPI by adding mpl.rcParams override after plot_SWOF import
"""
import json, pathlib

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "4. roughness_scale-cf-channel_eqs.ipynb"
nb = json.loads(NB.read_text())
cells = nb["cells"]

changes = 0

# ── 1. Fix Yen titles in _METHODS definition ──────────────────────────────
for ci, cell in enumerate(cells):
    src = "".join(cell["source"])
    if "_METHODS" in src and "'Yen R^{1/6}'" in src:
        new_lines = []
        for line in cell["source"]:
            line = line.replace("'Yen R^{1/6}'", r"'Yen $R^{1/6}$'")
            line = line.replace("'Yen R^{1/3}'", r"'Yen $R^{1/3}$'")
            new_lines.append(line)
        cell["source"] = new_lines
        changes += 1
        print(f"Fixed Yen titles in cell {ci}")
        break
else:
    print("WARNING: could not find _METHODS cell with Yen titles")

# ── 2. Increase display DPI ───────────────────────────────────────────────
DPI_LINE = "\n# Override low default DPI from plot_SWOF for sharper inline figures\nimport matplotlib as _mpl_dpi\n_mpl_dpi.rcParams['figure.dpi'] = 150\n"
for ci, cell in enumerate(cells):
    src = "".join(cell["source"])
    if "from source_functions_1p3 import" in src and "from plot_SWOF import" in src:
        # Check if DPI override already present
        if "figure.dpi" in src:
            print(f"DPI override already present in cell {ci}, skipping")
        else:
            cell["source"].append(DPI_LINE)
            changes += 1
            print(f"Added DPI=150 override in cell {ci}")
        break
else:
    print("WARNING: could not find plot_SWOF import cell")

if changes:
    NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
    print(f"Done — {changes} change(s) written.")
else:
    print("No changes needed.")
