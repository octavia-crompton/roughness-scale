#!/usr/bin/env python3
"""Patch cell 13: earthy palette + BrBG heatmap cmap."""
import json, pathlib

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"
nb = json.loads(NB.read_text("utf-8"))

TARGET = "Length-scale attribution of S_f decomposition terms"
cell = None
for c in nb["cells"]:
    src = "".join(c["source"])
    if TARGET in src:
        cell = c
        break
if cell is None:
    raise SystemExit(f"Could not find cell containing '{TARGET}'")

src = "".join(cell["source"])

# 1. Replace bold Dark2 palette with earthy tones
src = src.replace(
    "# -- Bold, saturated colour palette for all figures in this cell --\n"
    "_bold_colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']",
    "# -- Earthy colour palette for all figures in this cell --\n"
    "_earth_colors = ['#8c510a', '#bf812d', '#546223', '#35978f', '#8e6c3a', '#c45a28']",
)

# 2. Rename all remaining references
src = src.replace("_bold_colors", "_earth_colors")

# 3. Switch heatmap cmap from PuOr_r to BrBG (brown-teal, earthy diverging)
src = src.replace("cmap='PuOr_r'", "cmap='BrBG'")

cell["source"] = src.splitlines(keepends=True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print("OK – earthy palette applied")
