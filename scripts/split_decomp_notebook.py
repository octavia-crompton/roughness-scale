"""
Split roughness_scale-decomp.ipynb into:
  - roughness_scale-cf.ipynb      (setup + CF / composite-channel sections)
  - roughness_scale-decomp.ipynb  (setup + decomposition sections)

Split point: the first cell whose source contains '## Comparing effect-ratio predictions'.
"""

import json, re, copy
from pathlib import Path

NB_DIR = Path(__file__).resolve().parent.parent / "notebooks"
SRC    = NB_DIR / "roughness_scale-decomp.ipynb"

with open(SRC) as f:
    nb = json.load(f)

cells = nb["cells"]

# ── 1. Locate the split cell ──────────────────────────────────────────────────
SPLIT_TEXT = "## Comparing effect-ratio predictions"

split_idx = None
for i, cell in enumerate(cells):
    src = "".join(cell.get("source", []))
    if SPLIT_TEXT in src:
        split_idx = i
        break

if split_idx is None:
    raise RuntimeError(f"Could not find split marker: {SPLIT_TEXT!r}")

print(f"Split at cell index {split_idx}  (ID = {cells[split_idx].get('id', '?')})")

# ── 2. Identify setup cells (everything before the first ## Correction factor cell) ─
SETUP_TEXT = "## Correction factor"

setup_end = None
for i, cell in enumerate(cells):
    src = "".join(cell.get("source", []))
    if SETUP_TEXT in src:
        setup_end = i   # first non-setup cell
        break

if setup_end is None:
    raise RuntimeError(f"Could not find setup boundary marker: {SETUP_TEXT!r}")

print(f"Setup cells: 0..{setup_end - 1}  ({setup_end} cells)")
print(f"CF cells:    {setup_end}..{split_idx - 1}  ({split_idx - setup_end} cells)")
print(f"Decomp cells:{split_idx}..{len(cells) - 1}  ({len(cells) - split_idx} cells)")

setup_cells = cells[:setup_end]
cf_cells    = cells[setup_end:split_idx]
decomp_cells= cells[split_idx:]

# ── 3. Helper: patch notebook_name in setup cells ────────────────────────────
def patch_notebook_name(cells_list, new_name):
    """Replace the notebook_name string in the figure_registry.configure call."""
    patched = copy.deepcopy(cells_list)
    for cell in patched:
        src = "".join(cell.get("source", []))
        if "_fig_reg.configure" in src and "notebook_name=" in src:
            new_src = re.sub(
                r"notebook_name='[^']*'",
                f"notebook_name='{new_name}'",
                src,
            )
            # re-split into lines the same way Jupyter stores them
            lines = [l + "\n" for l in new_src.split("\n")]
            lines[-1] = lines[-1].rstrip("\n")   # no trailing newline on last line
            cell["source"] = lines
    return patched

setup_cf    = patch_notebook_name(setup_cells, "roughness_scale-cf.ipynb")
setup_decomp= patch_notebook_name(setup_cells, "roughness_scale-decomp.ipynb")

# ── 4. Build the two notebooks ───────────────────────────────────────────────
def make_nb(base, cells_list):
    new_nb = copy.deepcopy(base)
    new_nb["cells"] = cells_list
    # reset execution counts so re-run starts fresh
    for cell in new_nb["cells"]:
        if cell["cell_type"] == "code":
            cell["execution_count"] = None
            cell["outputs"] = []
    return new_nb

nb_cf    = make_nb(nb, setup_cf    + cf_cells)
nb_decomp= make_nb(nb, setup_decomp + decomp_cells)

# patch title markdown in CF notebook
for cell in nb_cf["cells"]:
    if cell["cell_type"] == "markdown":
        src = "".join(cell.get("source", []))
        if "Roughness scaling: analysis" in src:
            cell["source"] = [
                "## Roughness scaling: correction factor & composite channel equations\n",
                "\n",
                "_CF predictions, equivalent-roughness formulas (Lotter, Cox, Horton–Einstein, Felkel, …), and channel-equation comparisons._\n",
            ]
            break

# patch title markdown in decomp notebook
for cell in nb_decomp["cells"]:
    if cell["cell_type"] == "markdown":
        src = "".join(cell.get("source", []))
        if "Roughness scaling: analysis" in src:
            cell["source"] = [
                "## Roughness scaling: spatial decomposition\n",
                "\n",
                "_T0/T1/T2 variance decomposition of n_e/⟨n⟩, OLS and ML regression models, hybrid predictions, and pattern–forcing interactions._\n",
            ]
            break

# ── 5. Write outputs ──────────────────────────────────────────────────────────
CF_OUT    = NB_DIR / "roughness_scale-cf.ipynb"
DECOMP_OUT= NB_DIR / "roughness_scale-decomp.ipynb"

with open(CF_OUT, "w") as f:
    json.dump(nb_cf, f, indent=1)
print(f"Wrote {CF_OUT}  ({len(nb_cf['cells'])} cells)")

with open(DECOMP_OUT, "w") as f:
    json.dump(nb_decomp, f, indent=1)
print(f"Wrote {DECOMP_OUT}  ({len(nb_decomp['cells'])} cells)")

print("Done.")
