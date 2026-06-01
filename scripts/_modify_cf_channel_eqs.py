"""
Modify roughness_scale-cf-channel_eqs.ipynb:
  1. Replace the method-table markdown cell with expanded equations
  2. Replace the 2×2 plot cell → 1×3 (Leading term, T₀, T₁); drop Lotter/Felkel
  3. Fix fig5 → fig6 naming (from prior renumber that didn't persist)
"""
import json
from pathlib import Path

NB_PATH = Path("notebooks/roughness_scale-cf-channel_eqs.ipynb")

with open(NB_PATH) as f:
    nb = json.load(f)

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  Replace the markdown cell that has the method table
# ═══════════════════════════════════════════════════════════════════════════════
NEW_MARKDOWN = r"""### Perturbation expansion of the effect ratio

Starting from the friction slope $S_f = n^2 U^2 h^{-4/3}$, we write each field as a mean plus normalised fluctuation: $n = \bar{n}(1+r)$, $U = \bar{U}(1+\upsilon)$, $h = \bar{h}(1+\eta)$.

Spatially averaging and expanding to second order:

$$\frac{n_e}{\langle n\rangle} \approx 1 + \underbrace{\langle r^2\rangle + \langle \upsilon^2\rangle + 4\langle r\upsilon\rangle}_{T_0} \underbrace{- \tfrac{8}{3}\langle r\eta\rangle - \tfrac{8}{3}\langle \eta\upsilon\rangle}_{T_1 \text{ additions}} + \underbrace{\tfrac{14}{9}\langle \eta^2\rangle}_{T_2 \text{ addition}}$$

**Three approximation levels (general form):**

| Level | Formula |
|-------|---------|
| **Leading term** | $n_e/\langle n\rangle = \langle h\rangle^{2/3}\, S_0^{1/2}\; /\; \bigl(\langle U\rangle\,\langle n\rangle\bigr)$ |
| **$T_0$ ($r$, $\upsilon$ terms)** | $n_e/\langle n\rangle \approx 1 + \langle r^2\rangle + \langle \upsilon^2\rangle + 4\langle r\upsilon\rangle$ |
| **$T_1$ (adds $\eta$ cross-terms)** | $n_e/\langle n\rangle \approx T_0 - \tfrac{8}{3}\langle r\eta\rangle - \tfrac{8}{3}\langle \eta\upsilon\rangle$ |

**After CF substitution** — the kinematic-wave point solution gives $h_i \propto n_i^{3/5}$ and $U_i \propto n_i^{-3/5}$, so to first order $\eta \approx \tfrac{3}{5}r$ and $\upsilon \approx -\tfrac{3}{5}r$. All dispersive terms become multiples of $\langle r^2\rangle$:

| Level | CF-substituted form |
|-------|---------------------|
| **Leading term (CF)** | $n_e/\langle n\rangle = \langle h\rangle_{\rm CF}^{2/3}\, S_0^{1/2}\; /\; \bigl(\langle U\rangle_{\rm CF}\,\langle n\rangle\bigr)$ |
| **$T_0$ (CF)** | $n_e/\langle n\rangle \approx 1 - \tfrac{26}{25}\langle r^2\rangle$ |
| **$T_1$ (CF)** | $n_e/\langle n\rangle \approx 1 - \tfrac{42}{25}\langle r^2\rangle$ |

**CF coefficient derivation:**
- $\langle \upsilon^2\rangle = \tfrac{9}{25}\langle r^2\rangle$, $\;\;4\langle r\upsilon\rangle = -\tfrac{12}{5}\langle r^2\rangle$
  → $T_0$ coefficient: $1 + \tfrac{9}{25} - \tfrac{12}{5} = -\tfrac{26}{25}$
- $-\tfrac{8}{3}\langle r\eta\rangle = -\tfrac{8}{5}\langle r^2\rangle$, $\;\;-\tfrac{8}{3}\langle \eta\upsilon\rangle = +\tfrac{24}{25}\langle r^2\rangle$
  → $T_1$ adds $-\tfrac{16}{25}\langle r^2\rangle$, giving total $-\tfrac{42}{25}$
"""

md_found = False
for i, cell in enumerate(nb["cells"]):
    if cell.get("cell_type") != "markdown":
        continue
    src = "".join(cell.get("source", []))
    if "Equivalent roughness predictions" in src and "Leading term" in src:
        cell["source"] = NEW_MARKDOWN.splitlines(keepends=True)
        cell["outputs"] = []
        md_found = True
        print(f"  Replaced markdown cell at index {i}")
        break

if not md_found:
    print("  WARNING: markdown cell not found")

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  Replace the plot cell
# ═══════════════════════════════════════════════════════════════════════════════
NEW_PLOT = r"""# ── Predicted vs observed effect ratio — CF perturbation expansions ────────────
# Three approximation levels using only the roughness field (no simulation output).
# Markers are coloured by the variable that explains the most variance in the
# Leading-term (CF) residuals (pred − obs).

from stats import rmse as _rmse, r2 as _r2, best_residual_correlate

# ── observed column ─────────────────────────────────────────────────────────────
er_obs_col = 'effect_ratio'

# ── compute CF-based equivalent roughness predictions ──────────────────────────
# Leading term: direct back-calculation from CF-averaged h, U
summary['er_lead_CF'] = (
    summary['<h>_CF']**(2/3) * summary['So']**0.5
    / (summary['<U>_CF'] * summary['<n>'])
)

# T₀: keeps r, υ terms only  (υ ≈ -3/5 r under CF)
summary['er_T0_CF'] = (
    1
    + summary['<r2>']                        # ⟨r²⟩
    + (9/25) * summary['<r2>']               # ⟨υ²⟩ = (3/5)² ⟨r²⟩
    + 4 * (-3/5) * summary['<r2>']           # 4⟨rυ⟩ = 4·(-3/5)·⟨r²⟩
)

# T₁: adds η cross-terms  (η ≈ 3/5 r under CF)
summary['er_T1_CF'] = (
    1
    + summary['<r2>']                        # ⟨r²⟩
    + (9/25) * summary['<r2>']               # ⟨υ²⟩
    + 4 * (-3/5) * summary['<r2>']           # 4⟨rυ⟩
    - (8/3) * (3/5) * summary['<r2>']        # -8/3 ⟨rη⟩
    - (8/3) * (-9/25) * summary['<r2>']      # -8/3 ⟨ηυ⟩
)

# ── define 3-panel method list ──────────────────────────────────────────────────
methods = [
    ('er_lead_CF',   'Leading term (CF)'),
    ('er_T0_CF',     r'$T_0$: $r, \upsilon$ terms (CF)'),
    ('er_T1_CF',     r'$T_1$: adds $\eta$ cross-terms (CF)'),
]

# ── find variable that best explains Leading-term residuals ─────────────────────
_obs_arr    = summary[er_obs_col].to_numpy(float)
_lead_pred  = summary['er_lead_CF'].to_numpy(float)
_m_lead     = np.isfinite(_obs_arr) & np.isfinite(_lead_pred)
_lead_resid = _lead_pred[_m_lead] - _obs_arr[_m_lead]

_candidates = best_residual_correlate(
    summary, _lead_resid, mask=_m_lead,
    sim_cols=['fV', 'sigma', 'p', 'tr', 'l', 'aniso'],
)
COLOR_BY_VAR_4P   = _candidates[0][0] if _candidates else 'fV'
_pretty           = {'fV': r'$f_V$', 'sigma': r'$\sigma$', 'aniso': 'anisotropy',
                     'p': r'$p$ (mm/h)', 'tr': r'$t_r$ (min)', 'l': r'$L$ (m)'}
COLOR_BY_LABEL_4P = _pretty.get(COLOR_BY_VAR_4P, COLOR_BY_VAR_4P)
_cmap_name_4p     = VAR_CMAPS.get(COLOR_BY_VAR_4P, 'viridis')
_cmap_4p_base     = plt.get_cmap(_cmap_name_4p)
_cmap_4p          = mpl.colors.LinearSegmentedColormap.from_list(
    f'trunc_{_cmap_name_4p}', _cmap_4p_base(np.linspace(0.2, 0.95, 256)))

if _candidates:
    print(f"Coloring by: {COLOR_BY_VAR_4P}  (|r| with leading-term residual = {abs(_candidates[0][1]):.3f})")
    print("Top correlates:", [(c, round(r, 3)) for c, r in _candidates[:4]])
else:
    print(f"Coloring by: {COLOR_BY_VAR_4P}  (no valid correlates found; using default)")

# ── build discrete colour levels ────────────────────────────────────────────────
_color_vals_raw = pd.to_numeric(summary[COLOR_BY_VAR_4P], errors='coerce')
_unique_vals    = sorted(_color_vals_raw.dropna().unique())
MAX_LEVELS_4P   = 6
if len(_unique_vals) > MAX_LEVELS_4P:
    _binned      = pd.qcut(_color_vals_raw, q=MAX_LEVELS_4P, duplicates='drop')
    _bin_levels  = _binned.cat.categories.tolist()
    _color_vals  = _binned.cat.codes.astype(float)
    _color_vals[_binned.isna()] = np.nan
    _color_levels = list(range(len(_bin_levels)))
    _level_label  = lambda i: f'[{_bin_levels[i].left:.3g}, {_bin_levels[i].right:.3g}]'
else:
    _color_vals   = _color_vals_raw.to_numpy(float)
    _color_levels = _unique_vals
    _level_label  = lambda v: f'{v}'

_norm_4p    = mpl.colors.Normalize(vmin=min(_color_levels), vmax=max(_color_levels))
_palette_4p = {cv: _cmap_4p(_norm_4p(cv)) for cv in _color_levels}

# ── 1×3 plot ────────────────────────────────────────────────────────────────────
from matplotlib.ticker import MaxNLocator

fs_lab, fs_tck, fs_leg = 13, 12, 11
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
axes2_flat  = axes2.ravel()

for idx, (ax, (pred_col, title)) in enumerate(zip(axes2_flat, methods)):
    obs  = summary[er_obs_col].to_numpy(float)
    pred = (summary[pred_col].to_numpy(float)
            if pred_col in summary.columns else np.full(len(summary), np.nan))

    finite = np.isfinite(obs) & np.isfinite(pred)
    if finite.sum() == 0:
        ax.set_title(f"{title}\n(no data)", fontsize=fs_lab)
        continue

    lo = min(np.nanmin(obs[finite]), np.nanmin(pred[finite]))
    hi = max(np.nanmax(obs[finite]), np.nanmax(pred[finite]))
    ax.plot([lo, hi], [lo, hi], 'k-', lw=1, zorder=0)

    _rng = np.random.default_rng(42)
    for cv in _color_levels:
        sel = finite & (np.abs(_color_vals - cv) < 1e-9)
        if sel.sum() == 0:
            continue
        _jitter = _rng.normal(0, 0.005, sel.sum())
        ax.scatter(obs[sel], pred[sel] + _jitter,
                   color=_palette_4p[cv], marker='o',
                   s=25, alpha=0.80, edgecolors='none', zorder=3)

    rmse = _rmse(pred[finite], obs[finite])
    r2v  = _r2(pred[finite], obs[finite])
    ax.set_title(title, fontsize=fs_lab)
    ax.text(0.97, 0.05, f'RMSE={rmse:.3f}\n$R^2$={r2v:.2f}',
            transform=ax.transAxes, fontsize=fs_tck - 1,
            ha='right', va='bottom', zorder=5)
    ax.tick_params(labelsize=fs_tck)
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

    ax.set_xlabel(r"observed $n_e/\langle n\rangle$", fontsize=fs_lab)
    if idx == 0:
        ax.set_ylabel(r"predicted $n_e/\langle n\rangle$", fontsize=fs_lab)

    # ── collect per-method metrics for dynamic registry text ─────────────────
    if not hasattr(fig2, '_method_metrics'):
        fig2._method_metrics = {}
    fig2._method_metrics[pred_col] = dict(title=title, rmse=rmse, r2=r2v)

fig2.suptitle(r"Predicted vs observed $n_e/\langle n\rangle$"
              f"  —  coloured by {COLOR_BY_LABEL_4P}", fontsize=fs_lab + 1)
fig2.tight_layout(rect=[0, 0, 0.84, 1])

# ── legend ───────────────────────────────────────────────────────────────────────
leg_ax = fig2.add_axes([0.85, 0.12, 0.13, 0.74])
leg_ax.set_axis_off()
from matplotlib.lines import Line2D
_leg_handles = [Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=_palette_4p[cv], markersize=8, label=_level_label(cv))
                for cv in _color_levels]
leg_ax.legend(handles=_leg_handles, title=COLOR_BY_LABEL_4P, loc='center left',
              fontsize=fs_leg, title_fontsize=fs_leg, frameon=True, handlelength=1.2)

# ── save figure ──────────────────────────────────────────────────────────────────
_fig_dir, _, _ = _fig_dirs()
_fig6_name = 'fig6_pred_vs_obs_CF.png'
fig2.savefig(_os.path.join(_fig_dir, _fig6_name), dpi=300, bbox_inches='tight')

# ── build dynamic registry text from per-method metrics ──────────────────────
_mm = fig2._method_metrics
_interp_lines = []
for _pc, _ttl in methods:
    _m = _mm.get(_pc, {})
    _r2s = f"R²={_m['r2']:.2f}" if np.isfinite(_m.get('r2', np.nan)) else "R² undefined"
    _rms = f"RMSE={_m['rmse']:.3f}" if np.isfinite(_m.get('rmse', np.nan)) else ""
    _interp_lines.append(f"  - {_m.get('title', _pc)}: {_r2s}, {_rms}".rstrip(', '))

# Identify best method by R² (excluding NaN)
_valid = [(k, v) for k, v in _mm.items() if np.isfinite(v.get('r2', np.nan))]
if _valid:
    _best_key  = max(_valid, key=lambda x: x[1]['r2'])
    _worst_key = min(_valid, key=lambda x: x[1]['r2'])
    _best_title,  _best_r2  = _best_key[1]['title'],  _best_key[1]['r2']
    _worst_title, _worst_r2 = _worst_key[1]['title'], _worst_key[1]['r2']
    _concise_6 = (f"Predicted vs observed n_e/<n> for three CF perturbation expansions "
                  f"(Leading term, T0, T1), coloured by {COLOR_BY_VAR_4P}. "
                  f"{_best_title} performs best (R²={_best_r2:.2f}); "
                  f"{_worst_title} is weakest (R²={_worst_r2:.2f}).")
else:
    _concise_6 = f"Predicted vs observed n_e/<n> for three CF perturbation expansions, coloured by {COLOR_BY_VAR_4P}."

update_figure_registry(
    'fig6', _fig6_name,
    'Predicted vs observed n_e/<n> for three CF perturbation expansion levels '
    '(Leading term, T0 with r/υ terms, T1 adding η cross-terms).\n'
    f'Markers coloured by {COLOR_BY_VAR_4P} (highest |r| with leading-term residual).\n'
    'RMSE/R² annotated per panel.\n'
    '\n'
    'Per-method performance:\n' + '\n'.join(_interp_lines) + '\n'
    '\n'
    f'Residuals correlate most with {COLOR_BY_VAR_4P}; '
    f'dominant errors occur at extreme {COLOR_BY_VAR_4P} values.',
    concise=_concise_6)

plt.show()
"""

plot_found = False
for i, cell in enumerate(nb["cells"]):
    src = "".join(cell.get("source", []))
    if "er_lead_CF" in src and "methods = [" in src and "update_figure_registry" in src:
        cell["source"] = NEW_PLOT.splitlines(keepends=True)
        cell["outputs"] = []
        plot_found = True
        print(f"  Replaced plot cell at index {i}")
        break

if not plot_found:
    print("  WARNING: plot cell not found")

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  Write back
# ═══════════════════════════════════════════════════════════════════════════════
with open(NB_PATH, "w") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write("\n")

print("  Done — notebook written to disk")
