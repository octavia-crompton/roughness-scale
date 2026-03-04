# Roughness-Scale Project — AI Agent Instructions

> **Last updated:** 2026-03-02
> Read this file at the start of every new chat session.
> Update this file when new conventions, preferences, or project structure changes are established in chat.

---

## 1. Project Overview

This workspace studies roughness scaling in overland flow over heterogeneous (vegetated) hillslopes using SWOF simulations. The key quantity is the **effective roughness ratio** $n_e / \langle n \rangle$ — how the equivalent Manning's n compares to the spatial mean.

### Workspace layout

```
roughness-scale/
├── .github/copilot-instructions.md   ← THIS FILE
├── overleaf/                          ← git clone of https://git.overleaf.com/663fbe37ea3d10144699322a (push/pull to sync with Overleaf)
├── swof_code/                         ← Python modules (simulation I/O, plotting)
│   ├── source_functions_1p3.py        ← core simulation helpers
│   ├── plot_SWOF.py                   ← legacy plotting + `names` dict
│   ├── read_SWOF.py, write_SWOF.py
│   ├── plot_config.py, topo.py
│   └── ...
├── src/
│   ├── __init__.py
│   ├── labels.py                      ← SHARED labels, cmaps, font sizes (see §3)
│   └── figure_registry.py             ← SHARED figure registry helpers (see §3b)
├── notebooks/
│   ├── roughness_scale-analysis.ipynb ← original combined notebook (source for split; kept as reference)
│   ├── roughness_scale-pattern.ipynb  ← spatial pattern & storm characteristics (figs 1-3)
│   ├── roughness_scale-cf.ipynb       ← correction factor & composite channel equations (Lotter, Cox, etc.)
│   ├── roughness_scale-decomp.ipynb   ← spatial decomposition & prediction methods (T0/T1/T2, OLS, ML)
│   ├── roughness_scale-compute.ipynb  ← data processing & derived columns
│   └── archive/                       ← old/dated notebooks
├── figures/
│   └── <batch_case>/                  ← e.g. runaround_smooth/
│       ├── fig4_*.png, fig5_*.png     ← publication figures
│       ├── figure_registry.txt        ← full registry (auto-generated)
│       ├── figure_registry_concise.txt← concise registry (auto-generated)
│       └── scratch/                   ← exploratory figures (never registered)
```

### Key variables

| Python name | Meaning |
|---|---|
| `out_dir` | Path to the simulation batch (e.g. `…/Tests/runaround_smooth`) |
| `summary` | Main DataFrame — one row per simulation, loaded from `summary_slim.pkl` |
| `names` | Legacy label dict from `plot_SWOF.py` (imported via `from plot_SWOF import *`) |
| `rename` | Comprehensive LaTeX label dict (defined in `src/labels.py`) |
| `methods` | List of 4 tuples `(pred_col, title)` for CF prediction methods |

---

## 2. Notebook Conventions

### Active analysis notebooks

| Notebook | Purpose |
|---|---|
| `roughness_scale-pattern.ipynb` | Scatter/grid plots of σ, fV, anisotropy + storm characteristics vs $n_e/\langle n \rangle$ |
| `roughness_scale-cf.ipynb` | Correction factor, equivalent-roughness formulas (Lotter, Cox, Horton–Einstein, Felkel), and composite-channel comparisons |
| `roughness_scale-decomp.ipynb` | T0/T1/T2 variance decomposition of $n_e/\langle n \rangle$, OLS/ML regression models, hybrid predictions, pattern–forcing interactions |
| `roughness_scale-compute.ipynb` | Load raw SWOF output, compute `effect_ratio`, write `summary_slim.pkl` |

### Imports & sys.path

All analysis notebooks set `sys.path` to include `swof_code/` and use:
```python
from plot_SWOF import *        # brings `names` dict into scope
from source_functions_1p3 import *
```

For shared labels (see §3):
```python
import sys as _sys
_sys.path.insert(0, "/Users/octaviacrompton/Projects/roughness-scale/src")
from labels import (
    updates, rename, renameit,
    VEG_COLORS, VEG_LABELS,
    FS_LABEL, FS_TITLE, FS_TICK, FS_LEG,
    VAR_CMAPS,
    format_name as _format_name_raw,
)
# Wrap so callers don't need to pass `names` explicitly
def format_name(fld, updates=updates):
    return _format_name_raw(fld, names, updates=updates)
```

For the figure registry (see §4):
```python
import sys as _sys
_sys.path.insert(0, "/Users/octaviacrompton/Projects/roughness-scale/src")
import figure_registry as _fig_reg

# Re-expose convenience names
_fig_dirs              = _fig_reg._fig_dirs
update_figure_registry = _fig_reg.update_figure_registry

# Call once after out_dir is set — use the actual notebook filename:
_fig_reg.configure(out_dir, notebook_name='roughness_scale-pattern.ipynb')  # or -cf or -decomp
```

### Python environment

- Python: `/opt/homebrew/bin/python3` (base conda)
- No virtual env required; packages assumed installed (numpy, pandas, matplotlib, seaborn, statsmodels, sklearn, scipy)

### Cell structure preferences

- Keep one logical task per cell.
- Markdown cells before major sections.
- Use `summary.copy()` before in-place mutations (learned from `dislay_fit` bug).

### USE_HYDRO toggle

All three analysis notebooks support a `USE_HYDRO = False` flag. When `True`, `effect_ratio` and `effect` are overwritten with their hydrograph-based equivalents (`effect_ratio_hydro`, `effect_hydro`).

---

## 3. Shared Code — `src/labels.py`

All display labels, colour maps, and formatting helpers shared between notebooks live in `src/labels.py`. **Do not duplicate these in notebook cells.**

### Exports

| Name | Type | Purpose |
|---|---|---|
| `updates` | dict | Overrides for legacy `names` dict |
| `format_name(fld, names, updates)` | function | LaTeX-format a field name (takes `names` explicitly) |
| `rename` | dict | ~95 column → LaTeX label mappings |
| `renameit(name, mapping)` | function | Safe lookup: returns `name` if not in `mapping` |
| `VEG_COLORS` | dict | Canonical hex colours per veg_type |
| `VEG_LABELS` | dict | Display labels per veg_type |
| `FS_LABEL, FS_TITLE, FS_TICK, FS_LEG` | int | Global font sizes (12, 13, 12, 12) |
| `VAR_CMAPS` | dict | Per-variable canonical matplotlib colormaps |

### Adding new labels

When a new column appears in `summary`, add its LaTeX label to `src/labels.py` in the `rename` dict — not in the notebook.

### VAR_CMAPS colour assignments

Every simulation parameter has a fixed colourmap. Always use `VAR_CMAPS.get(var, 'viridis')` when colouring by a variable. Current assignments:

| Variable | Colourmap |
|---|---|
| `fV` | Greens |
| `sigma` | mako |
| `tr` | Blues |
| `p` | YlOrRd |
| `aniso` | coolwarm |
| `l` | Purples |
| `<n>` | GnBu |
| `_rain_cm` | YlOrBr |
| `p*tr` | PuBu |
| `r_ratio` | BrBG |
| `alpha_v` | Oranges |
| `alpha_b` | RdPu |
| `So` | copper |

---

## 3b. Shared Code — `src/figure_registry.py`

All figure-registry logic lives in `src/figure_registry.py`. **Do not duplicate these functions in notebook cells.**

### Exports

| Name | Purpose |
|---|---|
| `configure(out_dir, notebook_name)` | Initialise registry for a session — call once after `out_dir` is set |
| `_fig_dirs()` | Returns `(fig_dir, scratch_dir, registry_path)` |
| `update_figure_registry(fig_id, filename, description, concise)` | Write/update a figure entry in both registries |
| `_parse_registry(path)` | Parse existing registry → `{fig_id: entry}` |
| `_rewrite_registry(entries, path)` | Re-sort and write both registry files |

---

## 4. Figures & Figure Registry

### Two-tier figure system

1. **Publication figures** → saved to `figures/<batch>/`, named `fig<N>_<description>.png`
2. **Scratch/exploratory figures** → saved to `figures/<batch>/scratch/`, **never registered**

### Saving a publication figure

Every publication figure save must call `update_figure_registry()`:

```python
_fig_dir, _, _ = _fig_dirs()
_name = 'fig4_obs_vs_pred_re_6panel.png'
fig.savefig(_os.path.join(_fig_dir, _name), dpi=300, bbox_inches='tight')
update_figure_registry(
    'fig4', _name,
    description='Full multi-line description with interpretation...',
    concise='Two-sentence human summary. Key takeaway.')
```

### Saving a scratch figure

```python
_, _scratch, _ = _fig_dirs()
fig.savefig(_os.path.join(_scratch, 'name.png'), dpi=200, bbox_inches='tight')
# NO registry call
```

### update_figure_registry() behaviour

- Accepts: `fig_id`, `filename`, `description` (full), `concise` (2 sentences: what + interpretation)
- Automatically re-sorts entries: main figures (fig1 → figN) first, SI figures at the end
- Writes **two files** on every call:
  - `figure_registry.txt` — full detailed registry with metadata
  - `figure_registry_concise.txt` — short humanized version
- Both files include: creation date, source notebook name, figure save directory
- Each entry records: filename, update timestamp, source notebook, save directory

### Registry entry format (full)

```
### fig4 ###
File     : fig4_obs_vs_pred_re_6panel.png
Updated  : 2026-02-28 13:36
Notebook : notebooks/roughness_scale-analysis.ipynb
Saved in : figures/runaround_smooth
────────────────────────────────────────
<description with interpretation>
### end fig4 ###
```

### Registry entry format (concise)

```
Fig 4 — fig4_obs_vs_pred_re_6panel.png
  <Two sentences: description + interpretation.>
```

### Figure numbering

- Main figures: `fig1`, `fig2`, `fig3`, … (sequential, in paper order)
- SI figures: `SI1`, `SI2`, … (sorted after main figures)
- Numbering should be stable; don't renumber existing figures without explicit request

### Figure style rules

- Use `FS_LABEL`, `FS_TITLE`, `FS_TICK`, `FS_LEG` from `src/labels.py` for font sizes
- Use `VAR_CMAPS` when colouring by a simulation variable
- Use `VEG_COLORS` / `VEG_LABELS` for veg_type categories
- RMSE/R² annotations: `ax.text()` in bottom-right, not in title
- Legend style: prefer `Line2D` circle markers over `Patch` rectangles for scatter legends
- Truncate colormaps with darker lower bound (e.g., 0.4–0.95) to avoid too-light colours near white
- Add small y-jitter (`rng.normal(0, 0.005)`) when observed values cluster on discrete levels
- Use `renameit()` for axis labels whenever possible
- Multi-column grid plots: share y-axis between left and center columns (columns 0 and 1), keep right column independent

---

## 5. LaTeX Figure Compilations

Each batch folder in `figures/` gets its own LaTeX document that compiles all registered (non-scratch) figures with their concise captions.

### Workspace layout

```
latex/
├── .gitignore                         ← ignores .aux/.log/etc, keeps .tex/.pdf
├── build_figure_doc.py                ← parses concise registry → .tex → .pdf
├── watch_registry.sh                  ← fswatch auto-rebuild on registry change
└── runaround_smooth_figures.tex       ← generated (do not hand-edit)
```

### How it works

1. `build_figure_doc.py` reads `figures/<batch>/figure_registry_concise.txt`
2. Generates `latex/<batch>_figures.tex` with one `\figure` per registry entry
3. Compiles to PDF with `pdflatex` (if installed)
4. Captions come from the concise registry; empty captions fall back to filename

### Usage

```bash
# One-shot build (all batches, or specify one)
python3 latex/build_figure_doc.py
python3 latex/build_figure_doc.py runaround_smooth

# Auto-rebuild on registry changes (requires: brew install fswatch)
./latex/watch_registry.sh
```

### TeX installation

If `pdflatex` is not available, the `.tex` file is still written and the script exits cleanly. Install with:
```bash
brew install --cask mactex-no-gui
```

### Conventions

- **Do not hand-edit** the generated `.tex` files — they are overwritten on every build
- Figure labels match registry IDs: `\label{fig:fig4}`, `\label{fig:SI1}`, etc.
- One document per batch folder (`runaround_smooth_figures.tex`, etc.)
- The concise registry caption is used as `\caption{…}`; if empty, filename is shown
- New batch folders are auto-discovered when `build_figure_doc.py` runs with no args

---

## 6. Notebook File Management

- **Active notebooks** live in `notebooks/`.
- **Primary notebooks** for analysis are `roughness_scale-pattern.ipynb`, `roughness_scale-cf.ipynb`, and `roughness_scale-decomp.ipynb`.
- `roughness_scale-analysis.ipynb` is the original combined notebook kept as reference — do not make new changes to it.
- **Archived notebooks** go to `notebooks/archive/` with a date suffix if not already dated.
- Do not create new notebooks without being asked — modify existing ones.

---

## 7. Editing Preferences

- When editing `.ipynb` files, use `replace_string_in_file` with exact text matching (include 3–5 lines of context).
- For complex multi-line notebook cell edits that fail with text matching, fall back to a Python script using `json.load` / `json.dump`.
- Always use `summary.copy()` before any operation that modifies a subset in-place.
- Prefer `_private` names (leading underscore) for cell-local temporaries to avoid polluting the kernel namespace.

---

## 8. Updating This File

**This file should be updated whenever:**
- A new shared convention is established (colour, font, label, file layout)
- A new figure is added to the registry numbering scheme
- LaTeX manuscript conventions are established
- New shared code is added to `src/`
- File organization rules change

Add new sections or update existing ones. Keep the format consistent.
