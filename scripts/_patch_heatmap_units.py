#!/usr/bin/env python3
"""Patch heatmap: fractions instead of percents, auto-scale colorbar."""
import json, pathlib

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"
nb = json.loads(NB.read_text())

TARGET = "Length-scale attribution of S_f decomposition terms"

for cell in nb["cells"]:
    if cell["cell_type"] != "code":
        continue
    src = "".join(cell["source"])
    if TARGET not in src:
        continue

    # 1. Change annotation from percent to 2-decimal fraction
    src = src.replace("{_v:+.1%}", "{_v:+.2f}")

    # 2. Replace hardcoded vmin/vmax with data-driven symmetric range
    old_imshow = (
        "_im = ax_attr.imshow(_heat_data.values, cmap='RdBu_r', vmin=-0.6, vmax=1.0,\n"
        "                     aspect='auto', interpolation='nearest')"
    )
    new_imshow = (
        "_vabs = max(abs(_heat_data.values.min()), abs(_heat_data.values.max()))\n"
        "_vabs = np.ceil(_vabs * 10) / 10  # round up to nearest 0.1\n"
        "_im = ax_attr.imshow(_heat_data.values, cmap='RdBu_r', vmin=-_vabs, vmax=_vabs,\n"
        "                     aspect='auto', interpolation='nearest')"
    )
    if old_imshow not in src:
        print("ERROR: could not find imshow block")
        exit(1)
    src = src.replace(old_imshow, new_imshow)

    cell["source"] = src.splitlines(keepends=True)
    print("Patched")
    break

NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
print("Done")
