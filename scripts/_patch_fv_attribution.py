#!/usr/bin/env python3
"""Insert fV-grouped version of attribution + slope plots after the aniso-grouped cell."""
import json, pathlib

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"
nb = json.loads(NB.read_text("utf-8"))

# Find the aniso attribution cell (contains "Length-scale attribution of S_f")
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
# Length-scale attribution — grouped by vegetation fraction (fV)
# ============================================================================
import statsmodels.api as sm

# Reuse _sub, _Tcols, _labels, _n_terms, _sigma_to_LV from previous cell

# ============================================================================
# 1. Covariance attribution — by fV
# ============================================================================
_fv_levels = sorted(_sub['fV'].dropna().unique())
_attrib_fv_rows = []

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
    _attrib_fv_rows.append([f'$f_V={_fv:g}$'] + _fracs)

_attrib_fv = pd.DataFrame(_attrib_fv_rows, columns=['fV'] + _labels)
_mean_attrib_fv = _attrib_fv[_labels].mean()

print("--- Covariance attribution (% of across-sigma variance in S) ---")
print("  grouped by vegetation fraction:\n")
for _lab, _val in zip(_labels, _mean_attrib_fv):
    print(f"  {_lab:40s}  {_val:+6.1%}")
print(f"\n  {'Sum':40s}  {_mean_attrib_fv.sum():+6.1%}")

# -- Figure: horizontal grouped bar chart — grouped by term, bars = fV --
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
ax_attr_fv.set_title('Covariance attribution — by vegetation fraction', fontsize=FS_TITLE)
ax_attr_fv.legend(fontsize=FS_LEG - 1, frameon=True, loc='lower right')
ax_attr_fv.tick_params(labelsize=FS_TICK)
fig_attr_fv.tight_layout()
plt.show()

# ============================================================================
# 2. Slope decomposition — by fV
# ============================================================================
_slope_by_fv = {}
_se_by_fv    = {}
_total_by_fv = {}

for _fv in _fv_levels:
    _gf = _sub.loc[_sub['fV'] == _fv].copy()
    if _gf.shape[0] < 5:
        continue
    _fname = f'$f_V={_fv:g}$'
    # Use aniso dummies as controls
    _aniso_dum = pd.get_dummies(_gf['aniso'], prefix='aniso', drop_first=True, dtype=float)
    _Xf = sm.add_constant(pd.concat([_gf[['_LV']], _aniso_dum], axis=1).astype(float))

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

print("\n--- Slope decomposition: dT_k / dLV, by fV ---\n")
for _fname in _slope_fv_names:
    print(f"  {_fname}:")
    for _lab in _labels:
        _s = _slope_by_fv[_fname][_lab]
        _se = _se_by_fv[_fname][_lab]
        print(f"    {_lab:40s}  {_s:+.4f}  (+/-{_se:.4f})")
    print(f"    {'Total S slope':40s}  {_total_by_fv[_fname]:+.4f}\n")

# -- Figure: horizontal grouped bar chart of slopes by fV --
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
ax_slope_fv.set_title(r'Sensitivity of each term to $L_V$ by vegetation fraction',
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
print("\nSaved 2 fV-grouped figures to scratch/")
"""

new_cell = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": NEW_SOURCE.splitlines(keepends=True),
}

# Insert after the aniso attribution cell
nb["cells"].insert(idx + 1, new_cell)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print(f"OK – inserted fV-grouped attribution cell at position {idx + 1}")
