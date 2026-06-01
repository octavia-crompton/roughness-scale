"""Patch cell 12 (index 11) in decomp notebook: apply fV tick labels to 2×3 plot."""
import json

NB = "/Users/octaviacrompton/Projects/roughness-scale/notebooks/3. roughness_scale-decomp.ipynb"
with open(NB) as f:
    nb = json.load(f)

cell = nb["cells"][11]  # cell 12 (0-based index 11)
src = "".join(cell["source"])

# Verify this is the right cell
assert "Second-order" in src and "2, 3, figsize" in src, "Wrong cell!"

old = "    axes[i].tick_params(labelsize=FS_TICK + 1)\n\nfor ax in axes[3:]:"
new = (
    "    axes[i].tick_params(labelsize=FS_TICK + 1)\n\n"
    "_stride = max(1, len(_fv_levels) // 6)\n"
    "axes[0].set_xticks(_tick_pos[::_stride])\n"
    "axes[0].set_xticklabels(_tick_lbls[::_stride])\n"
    "for ax in axes[3:]:"
)

assert old in src, f"Search string not found in cell!"
src = src.replace(old, new, 1)

cell["source"] = src.splitlines(keepends=True)

with open(NB, "w") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write("\n")
print("Done — applied fV tick labels to cell 12")
