#!/usr/bin/env python3
"""Replace per-fV-level grouping with 3 fV bins in the fV attribution cell."""
import json, pathlib

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"
nb = json.loads(NB.read_text("utf-8"))

TARGET = "Length-scale attribution \u2014 grouped by vegetation fraction"
idx = None
for i, c in enumerate(nb["cells"]):
    src = "".join(c["source"])
    if TARGET in src:
        idx = i
        break
if idx is None:
    raise SystemExit(f"Could not find cell containing '{TARGET}'")

NEW_SOURCE = r"""# ============================================================================
# Length-scale attribution — grouped by vegetation fraction (fV bins)
# ============================================================================
import statsmodels.api as sm

# Reuse _sub, _Tcols, _labels, _n_terms, _sigma_to_LV from previous cell

# -- Bin fV into 3 groups for stable estimates --
_sub['_fV_bin'] = pd.qcut(_sub['fV'], q=3, duplicates='drop')
_fv_bin_sorted = sorted(_sub['_fV_bin'].dropna().unique(), key=lambda x: x.left)
_fv_bin_labels = [f'$f_V \\in ({iv.left:.2f},\\,{iv.right:.2f}]$' for iv in _fv_bin_sorted]
_fv_bin_map = dict(zip(_fv_bin_sorted, _fv_bin_labels))
_sub['_fV_label'] = _sub['_fV_bin'].map(_fv_bin_map)
_fv_groups_ordered = _fv_bin_labels  # low → high

print(f"fV bin sizes:")
for _bl in _fv_groups_ordered:
    _n = (_sub['_fV_label'] == _bl).sum()
    print(f"  {_bl}  n={_n}")

# ============================================================================
# 1. Covariance attribution — by fV bin
# ============================================================================
_attrib_fv_rows = []

for _fname in _fv_groups_ordered:
    _g = _sub.loc[_sub['_fV_label'] == _fname]
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
    _attrib_fv_rows.append([_fname] + _fracs)

_attrib_fv = pd.DataFrame(_attrib_fv_rows, columns=['fV'] + _labels)
_mean_attrib_fv = _attrib_fv[_labels].mean()

print("\n--- Covariance attribution (% of across-sigma variance in S) ---")
print("  grouped by vegetation fraction bin:\n")
for _lab, _val in zip(_labels, _mean_attrib_fv):
    print(f"  {_lab:40s}  {_val:+6.1%}")
print(f"\n  {'Sum':40s}  {_mean_attrib_fv.sum():+6.1%}")

# -- Figure: horizontal grouped bar chart — grouped by term, bars = fV bin --
_fv_groups = _attrib_fv['fV'].tolist()
_n_fv = len(_fv_groups)
_greens_cmap = plt.cm.Greens(np.linspace(0.3, 0.9, _n_fv))

_y_pos = np.arange(_n_terms)
_bar_h = 0.8 / _n_fv

fig_attr_fv, ax_attr_fv = plt.subplots(figsize=(9, 5.5))
for _fi, _fname in enumerate(_fv_groups):
    _offsets = _y_pos + (_fi - (_n_fv - 1) / 2) * _bar_h
    _vals = _attrib_fv.loc[_attrib_fv['fV'] == _fname, _labels].values.flatten()
    ax_attr_fv.barh(_offsets, _vals, height=_bar_h, label=_fname,
                    color=_greens_cmap[_fi], edgecolor='0.2', linewidth=0.6)

ax_attr_fv.axvline(0, color='k', lw=0.8)
ax_attr_fv.set_yticks(_y_pos)
ax_attr_fv.set_yticklabels(_labels, fontsize=FS_TICK)
ax_attr_fv.set_xlabel('fraction of Var$(S)$', fontsize=FS_LABEL)
ax_attr_fv.set_title('Covariance attribution — by vegetation fraction (binned)', fontsize=FS_TITLE)
ax_attr_fv.legend(fontsize=FS_LEG - 1, frameon=True, loc='lower right')
ax_attr_fv.tick_params(labelsize=FS_TICK)
fig_attr_fv.tight_layout()
plt.show()

# ============================================================================
# 2. Slope decomposition — by fV bin
# ============================================================================
_slope_by_fv = {}
_se_by_fv    = {}
_total_by_fv = {}

for _fname in _fv_groups_ordered:
    _gf = _sub.loc[_sub['_fV_label'] == _fname].copy()
    if _gf.shape[0] < 5:
        continue
    # Use aniso + storm dummies as controls
    _aniso_dum = pd.get_dummies(_gf['aniso'], prefix='aniso', drop_first=True, dtype=float)
    _storm_dum = pd.get_dummies(_gf[['p','tr']].astype(str).agg('_'.join, axis=1),
                                prefix='storm', drop_first=True, dtype=float)
    _Xf = sm.add_constant(pd.concat([_gf[['_LV']], _aniso_dum, _storm_dum], axis=1).astype(float))

    _sl, _se_sl = {}, {}
    for _k, _lab in enumerate(_labels):
        _y = _gf[_Tcols[_k]].to_numpy()
        _fit = sm.OLS(_y, _Xf).fit()
        _sl[_lab]    = _fit.params['_LV']
        _se_sl[_lab] = _fit.bse['_LV']
    _fitS = sm.OLS(_gf['_S'].to_numpy(), _Xf).fit()
    _slope_by_fv[_fname] = _sl
    _se_by_fv[_fname]    = _se_sl
    _total_by_fv[_fname] = _fitS.params['_LV']

_slope_fv_names = list(_slope_by_fv.keys())
_n_sf = len(_slope_fv_names)

print("\n--- Slope decomposition: dT_k / dLV, by fV bin ---\n")
for _fname in _slope_fv_names:
    print(f"  {_fname}:")
    for _lab in _labels:
        _s = _slope_by_fv[_fname][_lab]
        _se = _se_by_fv[_fname][_lab]
        print(f"    {_lab:40s}  {_s:+.4f}  (+/-{_se:.4f})")
    print(f"    {'Total S slope':40s}  {_total_by_fv[_fname]:+.4f}\n")

# -- Figure: horizontal grouped bar chart of slopes by fV bin --
_greens_slope = plt.cm.Greens(np.linspace(0.3, 0.9, _n_sf))
_y_pos = np.arange(_n_terms)
_bar_h = 0.8 / _n_sf

fig_slope_fv, ax_slope_fv = plt.subplots(figsize=(9, 5.5))
for _fi, _fname in enumerate(_slope_fv_names):
    _offsets = _y_pos + (_fi - (_n_sf - 1) / 2) * _bar_h
    _vals = [_slope_by_fv[_fname][_lab] for _lab in _labels]
    _ses  = [_se_by_fv[_fname][_lab] for _lab in _labels]
    ax_slope_fv.barh(_offsets, _vals, xerr=_ses, height=_bar_h, label=_fname,
                     color=_greens_slope[_fi], edgecolor='0.2', linewidth=0.6,
                     capsize=2, error_kw=dict(ecolor='0.3', lw=0.8))

ax_slope_fv.axvline(0, color='k', lw=0.8)
ax_slope_fv.set_yticks(_y_pos)
ax_slope_fv.set_yticklabels(_labels, fontsize=FS_TICK)
ax_slope_fv.set_xlabel(r'$\partial T_k\,/\,\partial L_V$  (m$^{-1}$)', fontsize=FS_LABEL)
ax_slope_fv.set_title(r'Sensitivity of each term to $L_V$ by vegetation fraction (binned)',
                      fontsize=FS_TITLE)
ax_slope_fv.legend(fontsize=FS_LEG - 1, frameon=True, loc='lower right')
ax_slope_fv.tick_params(labelsize=FS_TICK)
fig_slope_fv.tight_layout()
plt.show()

# -- Save to scratch --
_, _scratch, _ = _fig_dirs()
fig_attr_fv.savefig(_os.path.join(_scratch, 'lengthscale_cov_attribution_by_fV.png'),
                    dpi=200, bbox_inches='tight')
fig_slope_fv.savefig(_os.path.join(_scratch, 'lengthscale_slope_decomp_by_fV.png'),
                     dpi=200, bbox_inches='tight')
print("\nSaved 2 fV-binned figures to scratch/")
"""

nb["cells"][idx]["source"] = NEW_SOURCE.splitlines(keepends=True)
nb["cells"][idx]["outputs"] = []
nb["cells"][idx]["execution_count"] = None

NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print(f"OK – replaced cell {idx} with fV-binned version")
