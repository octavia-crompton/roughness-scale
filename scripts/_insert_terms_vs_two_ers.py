"""Insert new cell #VSC-bb4b8759 in 3. roughness_scale-decomp.ipynb:
6 decomposition terms vs two effect ratios (equiv & vol), coloured by ratio type.
"""
import json
from pathlib import Path

p = Path('notebooks/3. roughness_scale-decomp.ipynb')
nb = json.loads(p.read_text())

code = '''# ══════════════════════════════════════════════════════════════════════════════
# Each of the 6 second-order decomposition terms vs two observed effect ratios
# (equivalent: r_equiv5/<n>;  volume: n_vol/<n>), coloured by ratio type.
# ══════════════════════════════════════════════════════════════════════════════
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.stats import pearsonr as _pearsonr

# Make sure both effect-ratio columns exist on `summary`
if 'effect_ratio_vol' not in summary.columns:
    summary['effect_ratio_vol'] = summary['n_vol'] / summary['<n>']
if 'effect_ratio_equiv' not in summary.columns:
    summary['effect_ratio_equiv'] = summary['r_equiv5'] / summary['<n>']

_subset = summary.loc[subset.index]   # same filtered subset as Fig 3
_x_equiv = _subset['effect_ratio_equiv'].to_numpy(float)
_x_vol   = _subset['effect_ratio_vol'].to_numpy(float)

_er_specs = [
    (_x_equiv, r'equivalent  ($r_{\\rm eq5}/\\langle n\\rangle$)', '#1f77b4'),
    (_x_vol,   r'volume  ($n_{\\rm vol}/\\langle n\\rangle$)',     '#d62728'),
]

fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True)

for i, (ts, title) in enumerate(zip(_term_series, _term_titles)):
    ax = axes.flat[i]
    _y = ts.to_numpy(float)
    _txt_lines = []
    for _x, _lab, _col in _er_specs:
        ax.scatter(_x, _y, color=_col, s=35, alpha=0.7,
                   edgecolor='k', linewidth=0.3, label=_lab)
        _mask = _np_local.isfinite(_x) & _np_local.isfinite(_y)
        if _mask.sum() > 2:
            _r, _ = _pearsonr(_x[_mask], _y[_mask])
            _txt_lines.append(fr'$r_{{\\rm {_lab.split()[0]}}} = {_r:.2f}$')
    ax.axhline(0, color='grey', ls='--', lw=0.7, zorder=0)

    if _txt_lines:
        ax.text(0.97, 0.03, '\\n'.join(_txt_lines),
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=_FS_ALL - 2,
                bbox=dict(boxstyle='round,pad=0.25', fc='white',
                          ec='grey', alpha=0.8))

    ax.set_ylabel(title, fontsize=_FS_ALL)
    ax.tick_params(labelsize=_FS_ALL)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

for ax in axes[-1, :]:
    ax.set_xlabel(r'observed $n_e/\\langle n\\rangle$', fontsize=_FS_ALL)

_legend_handles = [
    Line2D([0], [0], marker='o', linestyle='None',
           markerfacecolor=_col, markeredgecolor='k',
           markeredgewidth=0.4, markersize=7, label=_lab)
    for _x, _lab, _col in _er_specs
]
fig.legend(
    handles=_legend_handles,
    title='effect-ratio definition',
    fontsize=_FS_ALL, title_fontsize=_FS_ALL,
    loc='center right',
    bbox_to_anchor=(1.0, 0.5),
    borderaxespad=0,
    frameon=True,
)

fig.suptitle(
    r"Second-order decomposition terms vs observed $n_e/\\langle n\\rangle$ "
    r"(equivalent vs volume definitions)",
    fontsize=_FS_ALL, y=1.0,
)
plt.tight_layout(rect=[0, 0, 0.84, 0.97])
plt.show()
'''

# locate target cell by id
target = None
for c in nb['cells']:
    if c.get('id') == 'f04f8aac':
        target = c
        break
assert target is not None and target['cell_type'] == 'code'
assert ''.join(target['source']).strip() == '', f"not empty: {target['source']!r}"

target['source'] = code.splitlines(keepends=True)
target['outputs'] = []
target['execution_count'] = None

p.write_text(json.dumps(nb, indent=1, ensure_ascii=False))
print(f"Updated cell bb4b8759 with {len(target['source'])} lines")
