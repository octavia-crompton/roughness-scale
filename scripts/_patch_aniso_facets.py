#!/usr/bin/env python3
"""Insert a new cell showing effect_ratio vs fV for each aniso level (3 subplots),
placed between the markdown cell and the attribution cell."""
import json, pathlib

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"
nb = json.loads(NB.read_text("utf-8"))

# Find the attribution cell
TARGET = "Length-scale attribution of S_f decomposition terms"
idx = None
for i, c in enumerate(nb["cells"]):
    src = "".join(c["source"])
    if TARGET in src:
        idx = i
        break
if idx is None:
    raise SystemExit(f"Could not find cell containing '{TARGET}'")

NEW_SOURCE = r"""# ============================================================================
# Length-scale sensitivity of n_e/<n> — faceted by anisotropy
# ============================================================================
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

_aniso_query = "hydro_err < 0.05 and p == 8 and tr == 60"
_sub_all = summary.query(_aniso_query).copy()

_aniso_vals = sorted(_sub_all['aniso'].dropna().unique())
_aniso_titles = {v: ('gradient-aligned' if v < 1 else
                     'isotropic'        if v == 1 else
                     'contour-aligned') for v in _aniso_vals}

# sigma colour scale (mako, same as fig 3)
_sigma_to_LV  = summary.groupby('sigma')['LV'].median().to_dict()
_sigma_levels = sorted(_sub_all['sigma'].dropna().unique())
_mako         = _cm.get_cmap('mako')
_sigma_norm   = _mcolors.Normalize(vmin=0, vmax=max(list(_sigma_levels) + [5.5]))
_sigma_palette = {s: _mako(_sigma_norm(s)) for s in _sigma_levels}
_n_levels      = len(_sigma_levels)

# fV positions
_fv_levels = sorted(_sub_all['fV'].dropna().unique())
_fv_to_pos = {v: float(i) for i, v in enumerate(_fv_levels)}
_fV_offset_step = 0.15
_sigma_offsets  = {s: (i - (_n_levels - 1) / 2) * _fV_offset_step
                   for i, s in enumerate(_sigma_levels)}

_tick_pos  = list(range(len(_fv_levels)))
_tick_lbls = [f'{v:g}' for v in _fv_levels]
_stride    = max(1, len(_fv_levels) // 5)
_tick_lbls = [lbl if i % _stride == 0 else '' for i, lbl in enumerate(_tick_lbls)]

fig_aniso, axs = plt.subplots(1, len(_aniso_vals), figsize=(5 * len(_aniso_vals), 4.5),
                               sharey=True)
if len(_aniso_vals) == 1:
    axs = [axs]

for _ai, _aval in enumerate(_aniso_vals):
    ax = axs[_ai]
    _g = _sub_all.loc[_sub_all['aniso'] == _aval]
    _fV_arr = _g['fV'].map(_fv_to_pos).to_numpy(float)
    _sigma_arr = _g['sigma'].to_numpy()
    _y_arr = _g['effect_ratio'].to_numpy(float)

    for s in _sigma_levels:
        _sel = _sigma_arr == s
        if _sel.any():
            ax.scatter(_fV_arr[_sel] + _sigma_offsets[s], _y_arr[_sel],
                       color=_sigma_palette[s], s=45, alpha=0.85, edgecolor='0.3',
                       linewidth=0.3)

    ax.axhline(1, color='grey', ls='--', lw=0.7, zorder=0)
    ax.set_title(_aniso_titles[_aval], fontsize=FS_TITLE)
    ax.set_xlabel(renameit('fV', rename), fontsize=FS_LABEL)
    ax.set_xticks(_tick_pos)
    ax.set_xticklabels(_tick_lbls, fontsize=FS_TICK)
    ax.tick_params(labelsize=FS_TICK)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

axs[0].set_ylabel(r'$n_e/\langle n\rangle$', fontsize=FS_LABEL)

def _lv_label(s):
    if s in _sigma_to_LV:
        return f'{_sigma_to_LV[s]:.0f} m'
    return fr"$\sigma={s}$"

_legend_handles = [
    Line2D([0], [0], marker='o', linestyle='None',
           markerfacecolor=_sigma_palette[s], markeredgecolor='k',
           markeredgewidth=0.4, markersize=7, label=_lv_label(s))
    for s in _sigma_levels
]
fig_aniso.legend(
    handles=_legend_handles,
    title=r'$L_V$',
    fontsize=FS_LEG, title_fontsize=FS_LEG,
    loc='center right',
    bbox_to_anchor=(1.0, 0.55),
    frameon=True,
)

fig_aniso.suptitle(
    r"Observed $n_e/\langle n\rangle$ vs $f_V$ by pattern orientation",
    fontsize=FS_TITLE, y=1.02,
)
plt.tight_layout()
plt.subplots_adjust(right=0.88)
plt.show()

# Save to scratch
_, _scratch, _ = _fig_dirs()
fig_aniso.savefig(_os.path.join(_scratch, 'effect_ratio_by_aniso.png'),
                  dpi=200, bbox_inches='tight')
print("Saved to scratch/")
"""

new_cell = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": NEW_SOURCE.splitlines(keepends=True),
}

nb["cells"].insert(idx, new_cell)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print(f"OK – inserted new cell at position {idx}")
