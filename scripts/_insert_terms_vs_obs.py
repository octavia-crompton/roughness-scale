"""Insert decomposition terms vs observed effect_ratio_hydro plot into cell 15
of notebooks/3. roughness_scale-decomp.ipynb (the empty cell after the Fig 3
3x3 decomposition cell)."""
import json
from pathlib import Path

p = Path('notebooks/3. roughness_scale-decomp.ipynb')
nb = json.loads(p.read_text())

code = '''# ══════════════════════════════════════════════════════════════════════════════
# Each of the 6 second-order decomposition terms vs observed effect_ratio_hydro
# (same subset / sigma colouring as Fig 3 above)
# ══════════════════════════════════════════════════════════════════════════════
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.stats import pearsonr as _pearsonr

_obs = subset['effect_ratio_hydro'].to_numpy(float)

fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True)

for i, (ts, title) in enumerate(zip(_term_series, _term_titles)):
    ax = axes.flat[i]
    _y = ts.to_numpy(float)
    for s in _sigma_levels:
        _sel = _sigma_codes == s
        if _sel.any():
            ax.scatter(_obs[_sel], _y[_sel],
                       color=_sigma_palette[s], s=40, alpha=0.8,
                       edgecolor='k', linewidth=0.3)
    ax.axhline(0, color='grey', ls='--', lw=0.7, zorder=0)

    _mask = _np_local.isfinite(_obs) & _np_local.isfinite(_y)
    if _mask.sum() > 2:
        _r, _ = _pearsonr(_obs[_mask], _y[_mask])
        ax.text(0.97, 0.03, fr'$r = {_r:.2f}$',
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=_FS_ALL,
                bbox=dict(boxstyle='round,pad=0.25', fc='white',
                          ec='grey', alpha=0.8))

    ax.set_ylabel(title, fontsize=_FS_ALL)
    ax.tick_params(labelsize=_FS_ALL)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

for ax in axes[-1, :]:
    ax.set_xlabel(r'observed $n_e/\\langle n\\rangle$ (hydrograph)',
                  fontsize=_FS_ALL)

_legend_handles = [
    Line2D([0], [0], marker='o', linestyle='None',
           markerfacecolor=_sigma_palette[s], markeredgecolor='k',
           markeredgewidth=0.4, markersize=7, label=_lv_label(s))
    for s in _sigma_levels
]
fig.legend(
    handles=_legend_handles,
    title=r'$L_V$ (patch scale)',
    fontsize=_FS_ALL, title_fontsize=_FS_ALL,
    loc='center right',
    bbox_to_anchor=(1.0, 0.5),
    borderaxespad=0,
    frameon=True,
)

fig.suptitle(
    r"Second-order decomposition terms vs observed $n_e/\\langle n\\rangle$ "
    r"(hydrograph-based)",
    fontsize=_FS_ALL, y=1.0,
)
plt.tight_layout(rect=[0, 0, 0.88, 0.97])
plt.show()
'''

target = nb['cells'][15]
assert target['cell_type'] == 'code'
assert ''.join(target['source']).strip() == '', f"cell 15 not empty: {target['source']!r}"

lines = code.splitlines(keepends=True)
target['source'] = lines
target['outputs'] = []
target['execution_count'] = None

p.write_text(json.dumps(nb, indent=1, ensure_ascii=False))
print(f"Updated cell 15 with {len(lines)} lines")
