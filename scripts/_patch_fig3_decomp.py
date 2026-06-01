"""Patch Figure 3 cell in roughness_scale-decomp.ipynb (Option A layout)."""
import json
from pathlib import Path

NB = Path('notebooks/3. roughness_scale-decomp.ipynb')

NEW_SRC = '''# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 — Dispersive stress terms and their relationship to spatial pattern
# ══════════════════════════════════════════════════════════════════════════════
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as _np_local

_term_titles = [
    r"$\\langle r^2 \\rangle$",
    r"$\\langle \\upsilon^2 \\rangle$",
    r"$4\\langle r\\,\\upsilon \\rangle$",
    r"$-\\frac{8}{3}\\langle r\\,\\eta \\rangle$",
    r"$-\\frac{8}{3}\\langle \\eta\\,\\upsilon \\rangle$",
    r"$\\frac{14}{9}\\langle \\eta^2 \\rangle$",
]
subset = summary.query("hydro_err < 0.05 and p == 5 and tr == 60 and aniso > 1")

_term_series = [subset['<r2>'], subset['<ups2>'], 4*subset['<r ups>'],
     -8/3*subset['<r eta>'], -8/3*subset['<eta ups>'], 14/9*subset['<eta2>']]

_top_titles = [
    r"observed $n_e/\\langle n\\rangle$",
    r"second-order ($T_2$)",
    r"Lotter",
]
_top_series = [subset['effect_ratio'], subset['effect_ratio_T2'], subset['effect_ratio_Q']]

# ── sigma colour scale — same mako palette as Fig 2 ─────────────────────────
_sigma_to_LV  = summary.groupby('sigma')['LV'].median().to_dict()
_sigma_levels = sorted(subset['sigma'].dropna().unique())
_mako         = _cm.get_cmap('mako')
_sigma_norm   = _mcolors.Normalize(vmin=0, vmax=max(list(_sigma_levels) + [5.5]))
_sigma_palette = {s: _mako(_sigma_norm(s)) for s in _sigma_levels}
_n_levels      = len(_sigma_levels)

# ── fV axis: uniform integer positions, small offset per sigma ───────────────
_fv_levels = sorted(subset['fV'].dropna().unique())
_fv_to_pos = {v: float(i) for i, v in enumerate(_fv_levels)}
_fV_vals   = subset['fV'].map(_fv_to_pos).to_numpy(float)
_tick_pos  = list(range(len(_fv_levels)))
_tick_lbls = [f'{v:g}' for v in _fv_levels]
_stride    = max(1, len(_fv_levels) // 5)
_tick_lbls = [lbl if i % _stride == 0 else '' for i, lbl in enumerate(_tick_lbls)]

_fV_offset_step = 0.15
_sigma_offsets  = {s: (i - (_n_levels - 1) / 2) * _fV_offset_step
                   for i, s in enumerate(_sigma_levels)}
_sigma_codes    = subset['sigma'].to_numpy()

# ── 3×3 figure: split into top (row 1) and bottom (rows 2–3) subfigures so
#    each block has its own suptitle and a clear vertical gap between them.
fig = plt.figure(figsize=(17, 12))
_subfigs = fig.subfigures(2, 1, height_ratios=[1, 2], hspace=0.12)
_axes_top = _subfigs[0].subplots(1, 3, sharex=True, sharey=True)
_axes_bot = _subfigs[1].subplots(2, 3, sharex=True, sharey=True)

axes = _np_local.empty((3, 3), dtype=object)
axes[0, :] = _axes_top
axes[1:, :] = _axes_bot

for col in range(1, 3):
    axes[0, col].tick_params(labelleft=False)
for row in range(1, 3):
    for col in range(1, 3):
        axes[row, col].tick_params(labelleft=False)

# Bumped axis-label font size (Option A: titles become y-labels)
_FS_LABEL_BIG = FS_LABEL + 4

def _scatter_panel(ax, ys, ylabel, hline_val=0):
    for s in _sigma_levels:
        _sel = _sigma_codes == s
        if _sel.any():
            ax.scatter(
                _fV_vals[_sel] + _sigma_offsets[s],
                ys.to_numpy(float)[_sel],
                color=_sigma_palette[s], s=40, alpha=0.8,
            )
    ax.axhline(hline_val, color='grey', ls='--', lw=0.7, zorder=0)
    ax.set_ylabel(ylabel, fontsize=_FS_LABEL_BIG)
    ax.tick_params(labelsize=FS_TICK)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

for ci, (ts, ylab) in enumerate(zip(_top_series, _top_titles)):
    _scatter_panel(axes[0, ci], ts, ylab, hline_val=1)

for i, (ts, ylab) in enumerate(zip(_term_series, _term_titles)):
    row, col = divmod(i, 3)
    _scatter_panel(axes[row + 1, col], ts, ylab, hline_val=0)

for ax in list(axes[0, :]) + list(axes[2, :]):
    ax.set_xlabel(renameit('fV', rename), fontsize=_FS_LABEL_BIG)
for ax in axes.ravel():
    ax.set_xticks(_tick_pos)
    ax.set_xticklabels(_tick_lbls, rotation=0, ha='right', fontsize=FS_TICK)

# ── Shared legend — same style as Fig 2 ─────────────────────────────────────
def _lv_label(s):
    if s in _sigma_to_LV:
        return f'{_sigma_to_LV[s]:.0f} m'
    return fr"$\\sigma={s}$"

_legend_handles = [
    Line2D([0], [0], marker='o', linestyle='None',
           markerfacecolor=_sigma_palette[s], markeredgecolor='k',
           markeredgewidth=0.4, markersize=7, label=_lv_label(s))
    for s in _sigma_levels
]
fig.legend(
    handles=_legend_handles,
    title=r'$L_V$ (patch scale)',
    fontsize=FS_LEG, title_fontsize=FS_LEG,
    loc='center right',
    bbox_to_anchor=(1.0, 0.8),
    borderaxespad=0,
    frameon=True,
)

_subfigs[0].suptitle(
    r"Effect ratio  $n_e/\\langle n\\rangle$:  observed vs. predictions",
    fontsize=FS_TITLE + 2, y=1.02,
)
_subfigs[1].suptitle(
    r"$S_f$ second-order decomposition terms"
    r"  ($r = n'/\\bar{n}$,  $\\upsilon = U'/\\bar{U}$,  $\\eta = h'/\\bar{h}$)",
    fontsize=FS_TITLE + 2, y=1.0,
)
plt.subplots_adjust(right=0.88)
plt.show()


# ── Save figure 3 ────────────────────────────────────────────────────────────
_fig_dir, _, _ = _fig_dirs()
_name = 'fig3_Sf_decomp_3x3.png'
fig.savefig(_os.path.join(_fig_dir, _name), dpi=300, bbox_inches='tight')
update_figure_registry(
    'fig3', _name,
    description=(
        r"3×3 panel figure showing the $S_f$ spatial decomposition at fixed "
        r"storm ($p=8$ mm/hr, $t_r=60$ min, anisotropic patterns only). "
        r"Row 1: observed $n_e/\\langle n\\rangle$, second-order prediction ($T_2$), "
        r"and Lotter prediction vs $f_V$, coloured by patch lengthscale $L_V$. "
        r"Rows 2–3: the six second-order fluctuation terms "
        r"($\\langle r^2\\rangle$, $\\langle\\upsilon^2\\rangle$, "
        r"$4\\langle r\\upsilon\\rangle$, $-\\frac{8}{3}\\langle r\\eta\\rangle$, "
        r"$-\\frac{8}{3}\\langle\\eta\\upsilon\\rangle$, "
        r"$\\frac{14}{9}\\langle\\eta^2\\rangle$) vs $f_V$. "
        r"Top row and bottom block have independent suptitles and a vertical "
        r"separator; math expressions are placed on y-axes (not titles) so "
        r"they read as variables rather than panel headings."
    ),
    concise=(
        r"$S_f$ decomposition figure (3×3): row 1 compares observed, $T_2$, and "
        r"Lotter effect ratios vs $f_V$; rows 2–3 show the six dispersive-flux "
        r"terms. Coloured by patch lengthscale $L_V$. Math expressions appear "
        r"as y-axis labels, with separate suptitles for the two blocks."
    ),
)
'''

nb = json.loads(NB.read_text())
patched = 0
for c in nb['cells']:
    if c['cell_type'] != 'code':
        continue
    src = ''.join(c.get('source', []))
    if 'Figure 3 — Dispersive stress terms' in src:
        # store as list of lines preserving trailing newlines (Jupyter convention)
        lines = NEW_SRC.splitlines(keepends=True)
        c['source'] = lines
        c['outputs'] = []
        c['execution_count'] = None
        patched += 1

assert patched == 1, f'expected 1 patch, got {patched}'
NB.write_text(json.dumps(nb, indent=1))
print(f'Patched {patched} cell(s).')
