#!/usr/bin/env python3
"""Insert cell: summary of effect_ratio sensitivity to anisotropy."""
import json, pathlib

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"
nb = json.loads(NB.read_text("utf-8"))

# Find the aniso facet cell (contains "effect_ratio_by_aniso")
TARGET = "effect_ratio_by_aniso"
idx = None
for i, c in enumerate(nb["cells"]):
    src = "".join(c["source"])
    if TARGET in src:
        idx = i
        break
if idx is None:
    raise SystemExit(f"Could not find cell containing '{TARGET}'")

NEW_SOURCE = r"""# ============================================================================
# Sensitivity of n_e/<n> to anisotropy — summary statistics & OLS
# ============================================================================
import statsmodels.api as sm

_sub_sens = summary.query("hydro_err < 0.05 and p == 8 and tr == 60").copy()

_aniso_map_s = {}
for _v in sorted(_sub_sens['aniso'].unique()):
    if _v < 1:
        _aniso_map_s[_v] = 'gradient-aligned'
    elif _v == 1:
        _aniso_map_s[_v] = 'isotropic'
    else:
        _aniso_map_s[_v] = 'contour-aligned'
_sub_sens['_aniso_label'] = _sub_sens['aniso'].map(_aniso_map_s)

# ── Summary table ────────────────────────────────────────────────────────────
_stats_rows = []
for _a in sorted(_sub_sens['aniso'].unique()):
    _g = _sub_sens.loc[_sub_sens['aniso'] == _a, 'effect_ratio']
    _stats_rows.append({
        'orientation': _aniso_map_s[_a],
        'aniso': _a,
        'n': len(_g),
        'mean': _g.mean(),
        'std': _g.std(),
        'min': _g.min(),
        'max': _g.max(),
        'range': _g.max() - _g.min(),
    })
_stats_df = pd.DataFrame(_stats_rows)
print("── Effect ratio summary by anisotropy ──\n")
print(_stats_df.to_string(index=False, float_format='{:.4f}'.format))

# ── OLS: effect_ratio ~ aniso (with fV + sigma controls) ────────────────────
_sub_sens['_LV'] = _sub_sens['sigma'].map(
    summary.groupby('sigma')['LV'].median().to_dict())
_fv_dum = pd.get_dummies(_sub_sens['fV'], prefix='fV', drop_first=True, dtype=float)
_sig_dum = pd.get_dummies(_sub_sens['sigma'], prefix='sig', drop_first=True, dtype=float)
_X_aniso = sm.add_constant(
    pd.concat([_sub_sens[['aniso']], _fv_dum, _sig_dum], axis=1))
_fit_aniso = sm.OLS(_sub_sens['effect_ratio'].to_numpy(), _X_aniso).fit()

print(f"\n── OLS: effect_ratio ~ aniso + fV + sigma FEs ──")
print(f"   aniso coef  = {_fit_aniso.params['aniso']:+.4f}  "
      f"(SE {_fit_aniso.bse['aniso']:.4f},  p = {_fit_aniso.pvalues['aniso']:.3g})")
print(f"   R²          = {_fit_aniso.rsquared:.3f}")
print(f"   n           = {int(_fit_aniso.nobs)}")

# ── Figure: box + strip plot of effect_ratio by orientation ──────────────────
_order = ['gradient-aligned', 'isotropic', 'contour-aligned']
_blues3 = plt.cm.Blues(np.linspace(0.35, 0.85, 3))

fig_sens, ax_sens = plt.subplots(figsize=(7, 4.5))
import seaborn as sns
_bp = sns.boxplot(data=_sub_sens, x='_aniso_label', y='effect_ratio',
                  order=_order, palette=_blues3, width=0.5, fliersize=0,
                  ax=ax_sens, boxprops=dict(edgecolor='0.2', linewidth=0.8),
                  medianprops=dict(color='k', lw=1.5))
sns.stripplot(data=_sub_sens, x='_aniso_label', y='effect_ratio',
              order=_order, color='0.25', size=3, alpha=0.5, jitter=0.15,
              ax=ax_sens)
ax_sens.axhline(1, color='grey', ls='--', lw=0.7)
ax_sens.set_xlabel('pattern orientation', fontsize=FS_LABEL)
ax_sens.set_ylabel(r'$n_e/\langle n\rangle$', fontsize=FS_LABEL)
ax_sens.set_title(
    f'Effect ratio by anisotropy  '
    r'(OLS $\beta_{aniso}$' + f' = {_fit_aniso.params["aniso"]:+.4f})',
    fontsize=FS_TITLE)
ax_sens.tick_params(labelsize=FS_TICK)
fig_sens.tight_layout()
plt.show()

# Save to scratch
_, _scratch, _ = _fig_dirs()
fig_sens.savefig(_os.path.join(_scratch, 'effect_ratio_aniso_sensitivity.png'),
                 dpi=200, bbox_inches='tight')
print("\nSaved to scratch/")
"""

new_cell = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": NEW_SOURCE.splitlines(keepends=True),
}

# Insert after the aniso facet cell
nb["cells"].insert(idx + 1, new_cell)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print(f"OK – inserted sensitivity summary cell at position {idx + 1}")
