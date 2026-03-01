"""
Shared display labels, colour maps, and formatting helpers for the
roughness-scale notebooks (analysis & compute).

Usage (in a notebook)
---------------------
    sys.path.insert(0, "/Users/octaviacrompton/Projects/roughness-scale/src")
    from labels import (
        updates, format_name, rename, renameit,
        VEG_COLORS, VEG_LABELS,
        FS_LABEL, FS_TITLE, FS_TICK, FS_LEG,
        VAR_CMAPS,
    )
"""

# ── format_name helper (legacy `names` dict) ─────────────────────────────────
# `updates`: overrides/extensions for the legacy `names` dict (from plot_SWOF).
# Used by format_name().
updates = {
    "effect_ratio" : "r_e/r_{avg}",
    'Delta_I_v'    : '\\Delta IF_V',
    'Ks_v'         : 'K_s',
    'K_avg'        : ' K_{avg}',
    'I_trend'      : 'dI/dt',
    'I_trend/Ks'   : 'K_s^{-1} dI/dt',
}


def format_name(fld, names, updates=updates):
    """Format a field name using the legacy *names* dict + *updates* overrides.

    Parameters
    ----------
    fld : str
        Column / variable name.
    names : dict
        The legacy ``names`` dict (typically ``from plot_SWOF import *``).
    updates : dict, optional
        Extra entries that override *names*.
    """
    names.update(updates)
    if fld in names:
        return '${0}$'.format(names[fld]).replace("\\\\", "\\").replace("D_", "\Delta \ ").replace("_cal", "_c")
    else:
        return '${0}$'.format(fld).replace("\\\\", "\\").replace("D_", "\Delta ").replace("_cal", "_c")


# ── rename: comprehensive LaTeX labels for seaborn / matplotlib axes ──────────
rename = {
    # ===== roughness / effect metrics =====
    "effect"             : r"$r_{\mathrm{avg}} - r_e$",
    "effect_ratio"       : r"$n_e/\langle n \rangle$",
    "effect_ratio_hydro" : r"$n_e^{\mathrm{hydro}}/\langle n \rangle$",
    "effect_ratio_geom"  : r"$r_e/r_{\mathrm{geom}}$",
    "r_equiv"            : r"$r_e$",
    "r_equiv5"           : r"$r_e$",
    "r_avg"              : r"$r_{\mathrm{avg}}$",
    "r_avg_off"          : r"$r_{\mathrm{avg}}$",
    "r_gmean"            : r"$r_{\mathrm{geom}}$",
    "r_ratio"            : r"$r_b/r_v$",
    "alpha_b"            : r"$r_b$",
    "alpha_v"            : r"$r_v$",

    # ===== controls / meta =====
    "fV"                 : r"$f_V$",
    "sigma"              : r"$\sigma$",
    "aniso"              : "anisotropy",
    "tr"                 : r"$t_R$",
    "l"                  : r"$L$",
    "seed"               : "seed",
    "radius"             : r"$r$",
    "Re"                 : r"$Re$",
    "Re_all"             : r"$Re$",

    # ===== hydraulics, means =====
    "<U>"                : r"$\langle U\rangle$",
    "<h>"                : r"$\langle h\rangle$",
    "<n>"                : r"$\langle n\rangle$",
    "<U>/<h>"            : r"$\langle U\rangle/\langle h\rangle$",

    # correction-factor predictions
    "<U>_CF"             : r"$\mathrm{CF}:\ \langle U\rangle$",
    "<h>_CF"             : r"$\mathrm{CF}:\ \langle h\rangle$",

    # cross-section integrals (area-weighted means)
    "<Ua>"               : r"$\langle U_a\rangle$",
    "<ha>"               : r"$\langle h_a\rangle$",
    "<na>"               : r"$\langle n_a\rangle$",

    # geometric means
    "n_gmean"            : r"$\langle n\rangle_{\mathrm{geom}}$",

    # ===== gradients / time-derivatives =====
    "d<U>/dt"            : r"$\partial_t \langle U\rangle$",
    "<dU/dx>"            : r"$\langle \partial_x U\rangle$",
    "<U dU/dx>"          : r"$\langle U\,\partial_x U\rangle$",
    "<U><dU/dx>"         : r"$\langle U\rangle\,\langle \partial_x U\rangle$",
    "<Up dUp/dx>"        : r"$\langle U'\,\partial_x U'\rangle$",
    "<dh/dx>"            : r"$\langle \partial_x h\rangle$",
    "<dhp/dx>"           : r"$\langle \partial_x h'\rangle$",

    # ===== variances & covariances =====
    "<Up2>"              : r"$\langle U'^2\rangle$",
    "<hp2>"              : r"$\langle h'^2\rangle$",
    "<np2>"              : r"$\langle n'^2\rangle$",
    "<np Up>"            : r"$\langle n' U'\rangle$",
    "<Up hp>"            : r"$\langle U' h'\rangle$",
    "<np hp>"            : r"$\langle n' h'\rangle$",

    # dimensionless variance ratios
    "c_var"              : r"$c=\langle U'^2\rangle/\langle U\rangle^2$",
    "k_nvar"             : r"$k=\langle n'^2\rangle/\langle n\rangle^2$",

    # ===== friction slope (direct and series pieces) =====
    "<Sf>"               : r"$\langle S_f\rangle$",
    "<Sf>_direct"        : r"$\langle S_f\rangle$",

    # T0 pieces
    "<Sf>_nbar2_Ubar2"   : r"$\langle h\rangle^{-4/3}\,\langle n\rangle^{2}\,\langle U\rangle^{2}$",
    "<Sf>_nbar2_Up2"     : r"$\langle h\rangle^{-4/3}\,\langle n\rangle^{2}\,\langle U'^2\rangle$",
    "<Sf>_Ubar2_np2"     : r"$\langle h\rangle^{-4/3}\,\langle U\rangle^{2}\,\langle n'^2\rangle$",
    "<Sf>_cross_nU"      : r"$4\,\langle h\rangle^{-4/3}\,\langle n\rangle\,\langle U\rangle\,\langle n'U'\rangle$",

    # T1 pieces (linear in h')
    "<Sf>_C_Uphp_lin"    : r"$-\frac{8}{3}\,\langle h\rangle^{-7/3}\,\langle n\rangle^{2}\,\langle U\rangle\,\langle U'h'\rangle$",
    "<Sf>_C_nphp_lin"    : r"$-\frac{8}{3}\,\langle h\rangle^{-7/3}\,\langle n\rangle\,\langle U\rangle^{2}\,\langle n'h'\rangle$",

    # T2 piece (h'^2)
    "<Sf>_h2"            : r"$\frac{14}{9}\,\langle h\rangle^{-10/3}\,\langle n\rangle^{2}\,\langle U\rangle^{2}\,\langle h'^2\rangle$",

    # totals / diagnostics
    "<Sf>_T0"            : r"$T_0$",
    "<Sf>_T1"            : r"$T_1$",
    "<Sf>_T2"            : r"$T_2$",
    "Tsum"               : r"$T_{0+1+2}$",

    # hybrid / series diagnostics (optional but included for completeness)
    "<Sf>_series2"       : r"$\mathrm{series}_2$",
    "<Sf>_hybrid_quad"   : r"$\mathrm{hybrid}_{\mathrm{quad}}$",
    "series2_err"        : r"$\Delta\,\mathrm{series}_2$",
    "hybrid_quad_err"    : r"$\Delta\,\mathrm{hybrid}_{\mathrm{quad}}$",
    "hybrid_quad_frac_series": r"$\mathrm{frac(series)}$",
    "<Sf>_direct_wet"    : r"$\langle S_f\rangle\ \mathrm{(wet)}$",
    "<Sf>_hyb_base"      : r"$\mathrm{hyb: base}$",
    "<Sf>_hyb_Up2"       : r"$\mathrm{hyb:}\ \langle U'^2\rangle$",
    "<Sf>_hyb_np2"       : r"$\mathrm{hyb:}\ \langle n'^2\rangle$",
    "<Sf>_hyb_cross"     : r"$\mathrm{hyb:}\ \langle n'U'\rangle$",
    "<Sf>_hyb_num_sum"   : r"$\mathrm{hyb: sum}$",
    "<Sf>_hyb_num_gap"   : r"$\mathrm{hyb: gap}$",
    "<Sf>_T0_hyb"        : r"$T_0\ (\mathrm{hyb})$",

    # ===== background parameters (if you plot them) =====
    "So"                 : r"$S_0$",
    "Ks_v"               : r"$K_s$",
    "p"                  : r"$p$",
    "i"                  : r"$i$",
    "dx"                 : r"$\Delta x$",
    "dt"                 : r"$\Delta t$",

    # ===== misc stats on h' =====
    "eta_max"            : r"$\max|h'/\langle h\rangle|$",
    "eta_mean"           : r"$\langle|h'/\langle h\rangle|\rangle$",
    "eta_rms"            : r"$\sqrt{\langle (h'/\langle h\rangle)^2\rangle}$",
    "wet"                : "wet cells",
    "eta"                : "series cells",
}


def renameit(name: str, mapping=rename):
    """Safe lookup for pretty labels."""
    return mapping.get(name, name)


# ── Consistent veg-type colour / label conventions ────────────────────────────
VEG_COLORS = {
    'rand':   '#2166ac',   # blue  – random patches
    'v_band': '#33a02c',   # green – along-slope (v) bands
    'band':   '#2166ac',   # blue  – along-slope bands (same group as v_band)
    'blob':   '#d6604d',   # orange-red – blob patches
}
VEG_LABELS = {
    'rand':   'patchy',
    'v_band': 'v-bands',
    'band':   'contour',
    'blob':   'blob',
}

# ── Global font sizes ─────────────────────────────────────────────────────────
FS_LABEL = 12   # axis label fontsize
FS_TITLE = 13   # subplot title fontsize
FS_TICK  = 12   # tick-label fontsize
FS_LEG   = 12   # legend text fontsize

# ── Per-variable canonical colour maps ────────────────────────────────────────
# Use these whenever colouring by a simulation parameter so plots are consistent.
VAR_CMAPS = {
    'fV'       : 'Greens',    # vegetation fraction     → greens
    'sigma'    : 'mako',      # roughness contrast      → mako (blue-green-yellow)
    'tr'       : 'Blues',     # storm duration          → blues
    'p'        : 'YlOrRd',    # rainfall intensity      → yellow-orange-red
    'aniso'    : 'coolwarm',  # anisotropy              → diverging red-blue
    'l'        : 'Purples',   # domain length           → purples
    '<n>'      : 'GnBu',      # mean Manning n          → green-blue
    '_rain_cm' : 'YlOrBr',    # total rain depth        → yellow-orange-brown
    'p*tr'     : 'PuBu',      # rainfall depth (p·tr)   → blue-purple
    'r_ratio'  : 'BrBG',      # roughness ratio nb/nv   → diverging brown-blue-green
    'alpha_v'  : 'Oranges',   # veg Manning n           → oranges
    'alpha_b'  : 'RdPu',      # bare Manning n          → red-purple
    'So'       : 'copper',    # slope                   → copper
}
