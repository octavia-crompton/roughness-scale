#!/usr/bin/env python3
"""Rewrite slope decomp (plot 2) to group by aniso with blues colormap."""
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

OLD_SLOPE = """\
# ============================================================================
# 2. Slope decomposition -- dT_k / d(LV) with fV fixed effects
# ============================================================================
import statsmodels.api as sm

_sub['_LV'] = _sub['sigma'].map(_sigma_to_LV)
_fv_dummies = pd.get_dummies(_sub['fV'], prefix='fV', drop_first=True, dtype=float)
_X_base = sm.add_constant(pd.concat([_sub[['_LV']], _fv_dummies], axis=1))

_slopes = {}
_slopes_se = {}
for _k, _lab in enumerate(_labels):
    _y = _sub[_Tcols[_k]].to_numpy()
    _fit = sm.OLS(_y, _X_base).fit()
    _slopes[_lab] = _fit.params['_LV']
    _slopes_se[_lab] = _fit.bse['_LV']

# Also fit total S
_fit_S = sm.OLS(_sub['_S'].to_numpy(), _X_base).fit()
_slope_S = _fit_S.params['_LV']

print("\\n--- Slope decomposition: dT_k / dLV (fV fixed effects) ---\\n")
for _lab in _labels:
    _s = _slopes[_lab]
    _se = _slopes_se[_lab]
    print(f"  {_lab:40s}  {_s:+.4f}  (+/-{_se:.4f})")
print(f"\\n  {'Sum of term slopes':40s}  {sum(_slopes.values()):+.4f}")
print(f"  {'Total S slope':40s}  {_slope_S:+.4f}")

# -- Figure: horizontal bar chart of slopes --
fig_slope, ax_slope = plt.subplots(figsize=(8, 4))
_sy = np.arange(_n_terms)
_slope_vals = [_slopes[_lab] for _lab in _labels]
_slope_ses  = [_slopes_se[_lab] for _lab in _labels]
ax_slope.barh(_sy, _slope_vals, xerr=_slope_ses, height=0.6,
              color=_earth_colors[:_n_terms], edgecolor='0.2', linewidth=0.6,
              capsize=3, error_kw=dict(ecolor='0.3', lw=1))
ax_slope.axvline(0, color='k', lw=0.8)
ax_slope.axvline(_slope_S, color='0.4', ls='--', lw=1,
                 label=rf'total $\\partial S/\\partial L_V = {_slope_S:.4f}$')
ax_slope.set_yticks(_sy)
ax_slope.set_yticklabels(_labels, fontsize=FS_TICK)
ax_slope.set_xlabel(r'$\\partial T_k\\,/\\,\\partial L_V$  (m$^{-1}$)', fontsize=FS_LABEL)
ax_slope.set_title(r'Sensitivity of each term to patch scale $L_V$ (fV fixed effects)',
                   fontsize=FS_TITLE)
ax_slope.legend(fontsize=FS_LEG, frameon=True, loc='lower right')
ax_slope.tick_params(labelsize=FS_TICK)
fig_slope.tight_layout()
plt.show()"""

NEW_SLOPE = """\
# ============================================================================
# 2. Slope decomposition -- dT_k / d(LV), grouped by anisotropy
# ============================================================================
import statsmodels.api as sm

_sub['_LV'] = _sub['sigma'].map(_sigma_to_LV)

# Fit per aniso level
_slope_by_aniso = {}   # {aniso_name: {label: slope}}
_se_by_aniso    = {}
_total_by_aniso = {}

for _a in _aniso_levels:
    _ga = _sub.loc[_sub['aniso'] == _a].copy()
    if _ga.shape[0] < 5:
        continue
    _aname = _aniso_map[_a]
    _fv_dum = pd.get_dummies(_ga['fV'], prefix='fV', drop_first=True, dtype=float)
    _Xa = sm.add_constant(pd.concat([_ga[['_LV']], _fv_dum], axis=1))

    _sl, _se_sl = {}, {}
    for _k, _lab in enumerate(_labels):
        _y = _ga[_Tcols[_k]].to_numpy()
        _fit = sm.OLS(_y, _Xa).fit()
        _sl[_lab]    = _fit.params['_LV']
        _se_sl[_lab] = _fit.bse['_LV']
    _fitS = sm.OLS(_ga['_S'].to_numpy(), _Xa).fit()
    _slope_by_aniso[_aname] = _sl
    _se_by_aniso[_aname]    = _se_sl
    _total_by_aniso[_aname] = _fitS.params['_LV']

_slope_aniso_names = list(_slope_by_aniso.keys())
_n_sa = len(_slope_aniso_names)

print("\\n--- Slope decomposition: dT_k / dLV, by anisotropy ---\\n")
for _aname in _slope_aniso_names:
    print(f"  {_aname}:")
    for _lab in _labels:
        _s = _slope_by_aniso[_aname][_lab]
        _se = _se_by_aniso[_aname][_lab]
        print(f"    {_lab:40s}  {_s:+.4f}  (+/-{_se:.4f})")
    print(f"    {'Total S slope':40s}  {_total_by_aniso[_aname]:+.4f}\\n")

# -- Figure: horizontal grouped bar chart of slopes by anisotropy --
_blues_slope = plt.cm.Blues(np.linspace(0.35, 0.85, _n_sa))
_y_pos = np.arange(_n_terms)
_bar_h = 0.8 / _n_sa

fig_slope, ax_slope = plt.subplots(figsize=(9, 5.5))
for _ai, _aname in enumerate(_slope_aniso_names):
    _offsets = _y_pos + (_ai - (_n_sa - 1) / 2) * _bar_h
    _vals = [_slope_by_aniso[_aname][_lab] for _lab in _labels]
    _ses  = [_se_by_aniso[_aname][_lab] for _lab in _labels]
    ax_slope.barh(_offsets, _vals, xerr=_ses, height=_bar_h, label=_aname,
                  color=_blues_slope[_ai], edgecolor='0.2', linewidth=0.6,
                  capsize=2, error_kw=dict(ecolor='0.3', lw=0.8))

ax_slope.axvline(0, color='k', lw=0.8)
ax_slope.set_yticks(_y_pos)
ax_slope.set_yticklabels(_labels, fontsize=FS_TICK)
ax_slope.set_xlabel(r'$\\partial T_k\\,/\\,\\partial L_V$  (m$^{-1}$)', fontsize=FS_LABEL)
ax_slope.set_title(r'Sensitivity of each term to $L_V$ by pattern orientation',
                   fontsize=FS_TITLE)
ax_slope.legend(fontsize=FS_LEG, frameon=True, loc='lower right')
ax_slope.tick_params(labelsize=FS_TICK)
fig_slope.tight_layout()
plt.show()"""

src = src.replace(OLD_SLOPE, NEW_SLOPE)

cell["source"] = src.splitlines(keepends=True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print("OK – slope decomp grouped by aniso with blues")
