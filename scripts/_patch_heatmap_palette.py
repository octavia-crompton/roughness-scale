#!/usr/bin/env python3
"""Patch cell 13 of roughness_scale-decomp.ipynb:
  - Broader subset (remove aniso > 1)
  - Group heatmap by anisotropy instead of fV
  - Use vivid PuOr_r diverging colormap (less pastel)
  - Black font for all heatmap annotations
"""
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

# -- 1. Broaden subset filter (remove aniso > 1) --
src = src.replace(
    'summary.query("hydro_err < 0.05 and p == 8 and tr == 60 and aniso >  1")',
    'summary.query("hydro_err < 0.05 and p == 8 and tr == 60")',
)
src = src.replace(
    'summary.query("hydro_err < 0.05 and p == 8 and tr == 60 and aniso > 1")',
    'summary.query("hydro_err < 0.05 and p == 8 and tr == 60")',
)

# -- 2. Replace fV grouping with aniso grouping --
OLD_GROUPING = """\
_fv_levels = sorted(_sub['fV'].dropna().unique())
_attrib_rows = []

for _fv in _fv_levels:
    _g = _sub.loc[_sub['fV'] == _fv]
    if _g.shape[0] < 3:
        continue
    _S = _g['_S'].to_numpy()
    _var_S = np.var(_S, ddof=1)
    if _var_S < 1e-15:
        continue
    _fracs = []
    for _tc in _Tcols:
        _cov = np.cov(_g[_tc].to_numpy(), _S)[0, 1]
        _fracs.append(_cov / _var_S)
    _attrib_rows.append([_fv] + _fracs)

_attrib = pd.DataFrame(_attrib_rows, columns=['fV'] + _labels)
_mean_attrib = _attrib[_labels].mean()

print("--- Covariance attribution (% of across-sigma variance in S) ---")
print("  averaged over fV levels:\\n")
for _lab, _val in zip(_labels, _mean_attrib):
    print(f"  {_lab:40s}  {_val:+6.1%}")
print(f"\\n  {'Sum':40s}  {_mean_attrib.sum():+6.1%}")

# -- Figure: stacked bar chart per fV + average --
_plot_df = _attrib.copy()
_avg_row = pd.DataFrame([['avg'] + _mean_attrib.tolist()], columns=['fV'] + _labels)
_plot_df = pd.concat([_plot_df, _avg_row], ignore_index=True)"""

NEW_GROUPING = """\
# -- Anisotropy labels for y-axis --
_aniso_map = {}
for _v in sorted(_sub['aniso'].unique()):
    if _v < 1:
        _aniso_map[_v] = 'gradient-aligned'
    elif _v == 1:
        _aniso_map[_v] = 'isotropic'
    else:
        _aniso_map[_v] = 'contour-aligned'

_aniso_levels = sorted(_sub['aniso'].dropna().unique())
_attrib_rows = []

for _a in _aniso_levels:
    _g = _sub.loc[_sub['aniso'] == _a]
    if _g.shape[0] < 3:
        continue
    _S = _g['_S'].to_numpy()
    _var_S = np.var(_S, ddof=1)
    if _var_S < 1e-15:
        continue
    _fracs = []
    for _tc in _Tcols:
        _cov = np.cov(_g[_tc].to_numpy(), _S)[0, 1]
        _fracs.append(_cov / _var_S)
    _attrib_rows.append([_aniso_map[_a]] + _fracs)

_attrib = pd.DataFrame(_attrib_rows, columns=['aniso'] + _labels)
_mean_attrib = _attrib[_labels].mean()

print("--- Covariance attribution (% of across-sigma variance in S) ---")
print("  averaged over anisotropy levels:\\n")
for _lab, _val in zip(_labels, _mean_attrib):
    print(f"  {_lab:40s}  {_val:+6.1%}")
print(f"\\n  {'Sum':40s}  {_mean_attrib.sum():+6.1%}")

# -- Figure: heatmap data per anisotropy + average --
_plot_df = _attrib.copy()
_avg_row = pd.DataFrame([['avg'] + _mean_attrib.tolist()], columns=['aniso'] + _labels)
_plot_df = pd.concat([_plot_df, _avg_row], ignore_index=True)"""

src = src.replace(OLD_GROUPING, NEW_GROUPING)

# -- 3. _soft_colors → keep for lollipop/profiles, but use vivid cmap for heatmap --
src = src.replace(
    "# -- Figure: annotated heatmap of covariance attribution --\n_heat_data = _plot_df.set_index('fV')[_labels].astype(float)",
    "# -- Figure: annotated heatmap of covariance attribution --\n_heat_data = _plot_df.set_index('aniso')[_labels].astype(float)",
)

# -- 4. Change heatmap colormap to PuOr_r (vivid, less pastel) --
src = src.replace("cmap='RdBu_r'", "cmap='PuOr_r'")

# -- 5. Black font for all annotations --
src = src.replace(
    "        _v = _heat_data.iloc[_i, _j]\n        _tc = 'w' if abs(_v) > 0.35 else '0.15'\n        ax_attr.text(_j, _i, f'{_v:+.2f}', ha='center', va='center',\n                     fontsize=FS_TICK, fontweight='bold', color=_tc)",
    "        _v = _heat_data.iloc[_i, _j]\n        ax_attr.text(_j, _i, f'{_v:+.2f}', ha='center', va='center',\n                     fontsize=FS_TICK, fontweight='bold', color='k')",
)

# -- 6. Y-labels: use aniso labels --
src = src.replace(
    "_ylabels_h = [f'$f_V={v:g}$' if isinstance(v, float) else v for v in _heat_data.index]",
    "_ylabels_h = list(_heat_data.index)",
)

# Write back
cell["source"] = src.splitlines(keepends=True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print("OK – patched cell successfully")
