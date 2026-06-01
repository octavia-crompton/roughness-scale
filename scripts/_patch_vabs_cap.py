"""Cap _vabs at 1.0 in the heatmap cell of notebook 3."""
import json, pathlib

nb_path = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"
nb = json.loads(nb_path.read_text())

old_line = "_vabs = np.ceil(_vabs * 10) / 10  # round up to nearest 0.1\n"
new_lines = [
    "_vabs = np.ceil(_vabs * 10) / 10  # round up to nearest 0.1\n",
    "_vabs = min(_vabs, 1.0)           # cap colour range so outliers don't wash out\n",
]

patched = False
for cell in nb["cells"]:
    src = cell.get("source", [])
    if old_line in src:
        idx = src.index(old_line)
        # Only patch if next line isn't already the cap
        if idx + 1 < len(src) and src[idx + 1].startswith("_vabs = min("):
            print("Already patched — skipping.")
            break
        cell["source"] = src[:idx] + new_lines + src[idx + 1:]
        patched = True
        print(f"Patched cell (inserted cap line after index {idx}).")
        break

if patched:
    nb_path.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
    print("Notebook written.")
else:
    print("Target line not found — no changes made.")
