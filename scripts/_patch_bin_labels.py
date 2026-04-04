"""Patch _bin_column in decomp notebook to use midpoint labels instead of ranges."""
import json, pathlib

nb_path = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

OLD = 'labels = [f"{iv.left:.2g}\u2013{iv.right:.2g}" for iv in cut.cat.categories]'
NEW = 'labels = [f"{(iv.left + iv.right) / 2:.2g}" for iv in cut.cat.categories]'

found = False
for cell in nb["cells"]:
    src = "".join(cell["source"])
    if OLD in src:
        new_src = src.replace(OLD, NEW)
        cell["source"] = new_src.splitlines(keepends=True)
        found = True
        break

if found:
    with open(nb_path, "w", encoding="utf-8") as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)
    print("OK – patched _bin_column labels to midpoints")
else:
    print("WARNING – old string not found; may already be patched")
