#!/usr/bin/env python3
"""Flip grouped bar: group by decomposition term (y-axis), bars = aniso levels. Horizontal."""
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

OLD = """\
# -- Figure: grouped bar chart of covariance attribution --
_groups = _plot_df['aniso'].tolist()
_n_groups = len(_groups)
_bar_width = 0.12
_x_pos = np.arange(_n_groups)

fig_attr, ax_attr = plt.subplots(figsize=(10, 5))
for _k, _lab in enumerate(_labels):
    _offsets = _x_pos + (_k - (_n_terms - 1) / 2) * _bar_width
    _vals = _plot_df[_lab].values
    ax_attr.bar(_offsets, _vals, width=_bar_width, label=_lab,
                color=_earth_colors[_k], edgecolor='0.2', linewidth=0.5)

ax_attr.axhline(0, color='k', lw=0.8)
# Separator line before the 'avg' group
ax_attr.axvline(_n_groups - 1.5, color='0.5', ls='--', lw=0.8)
ax_attr.set_xticks(_x_pos)
ax_attr.set_xticklabels(_groups, fontsize=FS_TICK)
ax_attr.set_ylabel('fraction of Var$(S)$', fontsize=FS_LABEL)
ax_attr.set_title('Covariance attribution of length-scale variance', fontsize=FS_TITLE)
ax_attr.legend(fontsize=FS_LEG - 1, frameon=True, ncol=2, loc='upper left')
ax_attr.tick_params(labelsize=FS_TICK)
fig_attr.tight_layout()
plt.show()"""

NEW = """\
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

src = src.replace(OLD, NEW)

cell["source"] = src.splitlines(keepends=True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print("OK – flipped to horizontal grouped bar by term")
