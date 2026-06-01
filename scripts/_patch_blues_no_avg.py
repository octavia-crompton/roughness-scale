#!/usr/bin/env python3
"""Drop 'avg' category from grouped bar, use blues colormap for aniso levels."""
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

# 1. Remove avg row construction
src = src.replace(
    "# -- Figure: heatmap data per anisotropy + average --\n"
    "_plot_df = _attrib.copy()\n"
    "_avg_row = pd.DataFrame([['avg'] + _mean_attrib.tolist()], columns=['aniso'] + _labels)\n"
    "_plot_df = pd.concat([_plot_df, _avg_row], ignore_index=True)",
    "# -- Figure: bar data per anisotropy (no average row) --\n"
    "_plot_df = _attrib.copy()",
)

# 2. Replace grouped bar section with blues colormap, no avg logic
OLD_BAR = """\
# -- Figure: horizontal grouped bar chart — grouped by term, bars = aniso level --
_aniso_groups = _plot_df['aniso'].tolist()   # e.g. ['gradient-aligned', 'isotropic', 'contour-aligned', 'avg']
_n_aniso = len(_aniso_groups)
_aniso_colors = ['#8c510a', '#bf812d', '#35978f', '#546223'][:_n_aniso]

_y_pos = np.arange(_n_terms)
_bar_h = 0.8 / _n_aniso

fig_attr, ax_attr = plt.subplots(figsize=(9, 5.5))
for _ai, _aname in enumerate(_aniso_groups):
    _offsets = _y_pos + (_ai - (_n_aniso - 1) / 2) * _bar_h
    _vals = _plot_df.loc[_plot_df['aniso'] == _aname, _labels].values.flatten()
    _lbl = _aname if _aname != 'avg' else 'average'
    _ec = '0.2' if _aname != 'avg' else '0.4'
    _hatch = '//' if _aname == 'avg' else None
    ax_attr.barh(_offsets, _vals, height=_bar_h, label=_lbl,
                 color=_aniso_colors[_ai], edgecolor=_ec, linewidth=0.6,
                 hatch=_hatch)

ax_attr.axvline(0, color='k', lw=0.8)
ax_attr.set_yticks(_y_pos)
ax_attr.set_yticklabels(_labels, fontsize=FS_TICK)
ax_attr.set_xlabel('fraction of Var$(S)$', fontsize=FS_LABEL)
ax_attr.set_title('Covariance attribution of length-scale variance', fontsize=FS_TITLE)
ax_attr.legend(fontsize=FS_LEG, frameon=True, loc='lower right')
ax_attr.tick_params(labelsize=FS_TICK)
fig_attr.tight_layout()
plt.show()"""

NEW_BAR = """\
# -- Figure: horizontal grouped bar chart — grouped by term, bars = aniso level --
_aniso_groups = _plot_df['aniso'].tolist()
_n_aniso = len(_aniso_groups)
_blues_cmap = plt.cm.Blues(np.linspace(0.35, 0.85, _n_aniso))

_y_pos = np.arange(_n_terms)
_bar_h = 0.8 / _n_aniso

fig_attr, ax_attr = plt.subplots(figsize=(9, 5.5))
for _ai, _aname in enumerate(_aniso_groups):
    _offsets = _y_pos + (_ai - (_n_aniso - 1) / 2) * _bar_h
    _vals = _plot_df.loc[_plot_df['aniso'] == _aname, _labels].values.flatten()
    ax_attr.barh(_offsets, _vals, height=_bar_h, label=_aname,
                 color=_blues_cmap[_ai], edgecolor='0.2', linewidth=0.6)

ax_attr.axvline(0, color='k', lw=0.8)
ax_attr.set_yticks(_y_pos)
ax_attr.set_yticklabels(_labels, fontsize=FS_TICK)
ax_attr.set_xlabel('fraction of Var$(S)$', fontsize=FS_LABEL)
ax_attr.set_title('Covariance attribution of length-scale variance', fontsize=FS_TITLE)
ax_attr.legend(fontsize=FS_LEG, frameon=True, loc='lower right')
ax_attr.tick_params(labelsize=FS_TICK)
fig_attr.tight_layout()
plt.show()"""

src = src.replace(OLD_BAR, NEW_BAR)

cell["source"] = src.splitlines(keepends=True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print("OK – dropped avg, switched to blues")
