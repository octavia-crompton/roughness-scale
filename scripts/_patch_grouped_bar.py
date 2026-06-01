#!/usr/bin/env python3
"""Replace the heatmap with a grouped bar plot (groups = anisotropy levels + avg)."""
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

OLD_HEATMAP = """\
# -- Figure: annotated heatmap of covariance attribution --
_heat_data = _plot_df.set_index('aniso')[_labels].astype(float)
_nrows_h, _ncols_h = _heat_data.shape

fig_attr, ax_attr = plt.subplots(figsize=(10, max(3.2, 0.55 * _nrows_h + 1.2)))
_vabs = max(abs(_heat_data.values.min()), abs(_heat_data.values.max()))
_vabs = np.ceil(_vabs * 10) / 10  # round up to nearest 0.1
_vabs = min(_vabs, 1.0)           # cap colour range so outliers don't wash out
_im = ax_attr.imshow(_heat_data.values, cmap='BrBG', vmin=-_vabs, vmax=_vabs,
                     aspect='auto', interpolation='nearest')
for _i in range(_nrows_h):
    for _j in range(_ncols_h):
        _v = _heat_data.iloc[_i, _j]
        ax_attr.text(_j, _i, f'{_v:+.2f}', ha='center', va='center',
                     fontsize=FS_TICK, fontweight='bold', color='k')
# separator before avg row
ax_attr.axhline(_nrows_h - 1.5, color='k', lw=1.2)
ax_attr.set_xticks(range(_ncols_h))
ax_attr.set_xticklabels(_labels, fontsize=FS_TICK - 1, rotation=25, ha='right')
_ylabels_h = list(_heat_data.index)
ax_attr.set_yticks(range(_nrows_h))
ax_attr.set_yticklabels(_ylabels_h, fontsize=FS_TICK)
_cb = plt.colorbar(_im, ax=ax_attr, shrink=0.8, pad=0.02)
_cb.set_label('fraction of Var$(S)$', fontsize=FS_LABEL)
_cb.ax.tick_params(labelsize=FS_TICK)
ax_attr.set_title('Covariance attribution of length-scale variance', fontsize=FS_TITLE)
fig_attr.tight_layout()
plt.show()"""

NEW_GROUPED_BAR = """\
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

src = src.replace(OLD_HEATMAP, NEW_GROUPED_BAR)

cell["source"] = src.splitlines(keepends=True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print("OK – replaced heatmap with grouped bar plot")
