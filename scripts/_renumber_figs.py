#!/usr/bin/env python3
"""
Renumber figures: insert new fig3 in pattern notebook, bump fig3→fig4, fig4→fig5, fig5→fig6.

Steps:
  1. Rename PNG files in figures/runaround_smooth/
  2. Renumber cells in decomp notebook (fig3 → fig4)
  3. Renumber cells in cf notebook (fig4 → fig5, fig5 → fig6)
  4. Renumber cells in cf-channel_eqs notebook (fig5 → fig6)
  5. Insert new fig3 save cell in pattern notebook after the LV figure cell
  6. Rewrite figure_registry.txt and figure_registry_concise.txt
"""
import json, os, re, shutil
from pathlib import Path
from datetime import datetime

ROOT = Path("/Users/octaviacrompton/Projects/roughness-scale")
FIG_DIR = ROOT / "figures" / "runaround_smooth"
NB_DIR = ROOT / "notebooks"

# ── 1. Rename PNG files (must go high→low to avoid collisions) ────────────
renames = [
    ("fig5_pred_vs_obs_CF.png", "fig6_pred_vs_obs_CF.png"),
    ("fig4_obs_vs_pred_re_6panel.png", "fig5_obs_vs_pred_re_6panel.png"),
    ("fig3_Sf_decomp_3x3.png", "fig4_Sf_decomp_3x3.png"),
]
for old, new in renames:
    src = FIG_DIR / old
    dst = FIG_DIR / new
    if src.exists():
        print(f"  rename: {old} → {new}")
        src.rename(dst)
    else:
        print(f"  SKIP (not found): {old}")

# ── helper: bulk string replace in a notebook cell's source list ──────────
def nb_replace(nb_path, replacements):
    """
    replacements: list of (old_str, new_str) pairs.
    Applies all replacements to every cell source line. Returns change count.
    """
    with open(nb_path, "r") as f:
        nb = json.load(f)
    changes = 0
    for cell in nb["cells"]:
        new_src = []
        for line in cell["source"]:
            orig = line
            for old, new in replacements:
                line = line.replace(old, new)
            if line != orig:
                changes += 1
            new_src.append(line)
        cell["source"] = new_src
    with open(nb_path, "w") as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
        f.write("\n")
    return changes

# ── 2. decomp notebook: fig3 → fig4 ──────────────────────────────────────
decomp = NB_DIR / "3. roughness_scale-decomp.ipynb"
n = nb_replace(decomp, [
    ("fig3_Sf_decomp_3x3.png", "fig4_Sf_decomp_3x3.png"),
    ("'fig3'", "'fig4'"),
])
print(f"  decomp: {n} line(s) changed")

# ── 3. cf notebook: fig4 → fig5, fig5 → fig6 (order matters!) ────────────
cf = NB_DIR / "roughness_scale-cf.ipynb"
# Do fig5→fig6 FIRST, then fig4→fig5 to avoid double-bumping
n = nb_replace(cf, [
    # fig5 → fig6
    ("_fig5_name", "_fig6_name"),
    ("fig5_pred_vs_obs_CF.png", "fig6_pred_vs_obs_CF.png"),
    ("'fig5'", "'fig6'"),
    # fig4 → fig5
    ("_fig4_name", "_fig5_name"),
    ("_fig4_metrics", "_fig5_metrics"),
    ("fig4_obs_vs_pred_re_6panel.png", "fig5_obs_vs_pred_re_6panel.png"),
    ("'fig4'", "'fig5'"),
    ("_valid4", "_valid5"),
    ("_m4d", "_m5d"),
    ("_r2_4", "_r2_5"),
    ("_rmse4", "_rmse5"),
])
print(f"  cf: {n} line(s) changed")

# ── 4. cf-channel_eqs notebook: fig5 → fig6 ──────────────────────────────
cf_ch = NB_DIR / "roughness_scale-cf-channel_eqs.ipynb"
n = nb_replace(cf_ch, [
    ("_fig5_name", "_fig6_name"),
    ("fig5_pred_vs_obs_CF.png", "fig6_pred_vs_obs_CF.png"),
    ("'fig5'", "'fig6'"),
])
print(f"  cf-channel_eqs: {n} line(s) changed")

# ── 5. Insert fig3 save cell in pattern notebook after the LV figure cell ─
pattern_nb = NB_DIR / "2. roughness_scale-pattern.ipynb"
with open(pattern_nb, "r") as f:
    nb = json.load(f)

# Find the cell whose id is VSC-8b78d522 (cell 22, the LV vs effect_ratio figure)
target_id = "VSC-8b78d522"
insert_idx = None
for i, cell in enumerate(nb["cells"]):
    if cell.get("id") == target_id:
        insert_idx = i + 1
        break

if insert_idx is None:
    print("  ERROR: could not find target cell in pattern notebook!")
else:
    save_cell = {
        "cell_type": "code",
        "metadata": {},
        "outputs": [],
        "source": [
            "# ── Save fig3 ────────────────────────────────────────────────────────\n",
            "_fig_dir, _, _ = _fig_dirs()\n",
            "_name = 'fig3_effect_ratio_vs_LV.png'\n",
            "fig.savefig(_os.path.join(_fig_dir, _name), dpi=300, bbox_inches='tight')\n",
            "update_figure_registry(\n",
            "    'fig3', _name,\n",
            "    description=(\n",
            "        '2×3 panel figure showing $n_e/\\\\langle n\\\\rangle$ vs patch lengthscale $L_V$ '\n",
            "        'for six vegetation fractions ($f_V = 0.1$ to $0.9$, excluding 0.05), '\n",
            "        'at a single storm ($p$, $t_r$). Points coloured and shaped by orientation '\n",
            "        '(gradient-aligned, isotropic, contour-aligned) with small horizontal offsets '\n",
            "        'for legibility. Shows that $n_e/\\\\langle n\\\\rangle$ generally increases with $L_V$, '\n",
            "        'approaching unity for the largest patches, and that the orientation effect '\n",
            "        'is strongest at intermediate $f_V$.'\n",
            "    ),\n",
            "    concise=(\n",
            "        'Effect ratio vs patch lengthscale $L_V$, panelled by $f_V$, coloured by orientation. '\n",
            "        'Ratio increases with patch size toward unity; orientation sensitivity peaks near $f_V \\\\approx 0.5$.'\n",
            "    ),\n",
            ")\n",
        ],
    }
    nb["cells"].insert(insert_idx, save_cell)
    with open(pattern_nb, "w") as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
        f.write("\n")
    print(f"  pattern: inserted fig3 save cell at position {insert_idx}")

# ── 6. Rewrite figure registries ─────────────────────────────────────────
now = datetime.now().strftime("%Y-%m-%d %H:%M")

# ---- Full registry ----
full_reg = f"""\
Figure registry  •  created {now}
Source notebook  :  notebooks/roughness_scale-pattern.ipynb
Figures saved in :  /Users/octaviacrompton/Projects/roughness-scale/figures/runaround_smooth
────────────────────────────────────────────────────────────────────────

### fig1 ###
File     : fig1_conceptual_sketch.png
Updated  : 2026-03-31 17:36
Notebook : notebooks/roughness_scale-compute.ipynb
Saved in : ../figures/runaround_smooth
────────────────────────────────────────────────────────────────────────
Conceptual figure: (a) 3-D view of a patchy hillslope with vegetation fraction $f_V=0.4$, $n_{{\\rm veg}}=0.20$, $n_{{\\rm bare}}=0.05$. (b) Stack of uniform-roughness simulations (library), colour-coded from high $n$ (dark green) to low $n$ (light green), with the best-match layer highlighted. (c) Outlet hydrograph: heterogeneous simulation (blue), uniform-roughness library (green shades), and matched $n_e$ (green dashed). The figure poses the central question: what equivalent $n$ reproduces the heterogeneous hydrograph?
### end fig1 ###

### fig2 ###
File     : fig2_effect_ratio_fV_2x2.png
Updated  : 2026-04-04 18:53
Notebook : notebooks/roughness_scale-pattern.ipynb
Saved in : ../figures/runaround_smooth
────────────────────────────────────────────────────────────────────────
2×2 panel figure. Top row: $n_e/\\langle n\\rangle$ vs vegetation fraction $f_V$ for (left) least intense and (right) most intense storm, coloured by patch lengthscale $L_V$ (mako palette). Dashed line at unity. Small horizontal offsets per $L_V$ level aid legibility. Bottom row: outlet hydrographs at $f_V\\approx 0.5$ for the same two storms, coloured by $L_V$. Shows that $n_e/\\langle n\\rangle < 1$ for all cases, with the minimum near $f_V\\approx 0.5$; ratio increases toward unity with larger patch scale; modest storm sensitivity.
### end fig2 ###

### fig3 ###
File     : fig3_effect_ratio_vs_LV.png
Updated  : {now}
Notebook : notebooks/2. roughness_scale-pattern.ipynb
Saved in : ../figures/runaround_smooth
────────────────────────────────────────────────────────────────────────
2×3 panel figure showing $n_e/\\langle n\\rangle$ vs patch lengthscale $L_V$ for six vegetation fractions ($f_V = 0.1$ to $0.9$, excluding 0.05), at a single storm ($p$, $t_r$). Points coloured and shaped by orientation (gradient-aligned, isotropic, contour-aligned) with small horizontal offsets for legibility. Shows that $n_e/\\langle n\\rangle$ generally increases with $L_V$, approaching unity for the largest patches, and that the orientation effect is strongest at intermediate $f_V$.
### end fig3 ###

### fig4 ###
File     : fig4_Sf_decomp_3x3.png
Updated  : 2026-04-03 18:40
Notebook : notebooks/roughness_scale-decomp.ipynb
Saved in : ../figures/runaround_smooth
────────────────────────────────────────────────────────────────────────
3×3 panel figure showing the $S_f$ spatial decomposition at fixed storm ($p=8$ mm/hr, $t_r=60$ min, anisotropic patterns only). Row 1: observed $n_e/\\langle n\\rangle$, second-order prediction ($T_2$), and Lotter prediction vs $f_V$, coloured by patch lengthscale $L_V$. Rows 2–3: the six second-order fluctuation terms ($\\langle r^2\\rangle$, $\\langle\\upsilon^2\\rangle$, $4\\langle r\\upsilon\\rangle$, $-\\tfrac{{8}}{{3}}\\langle r\\eta\\rangle$, $-\\tfrac{{8}}{{3}}\\langle\\eta\\upsilon\\rangle$, $\\tfrac{{14}}{{9}}\\langle\\eta^2\\rangle$) vs $f_V$. Row 1 shares a y-axis; rows 2–3 share a common y-axis to allow direct comparison of term magnitudes.
### end fig4 ###

### fig5 ###
File     : fig5_obs_vs_pred_re_6panel.png
Updated  : 2026-03-12 11:41
Notebook : notebooks/roughness_scale-cf.ipynb
Saved in : ../figures/runaround_smooth
────────────────────────────────────────────────────────────────────────
Observed vs predicted r_e for six equivalent-roughness expressions (Horton-Einstein, Pavlovskii, Lotter, Cox, Shear-force, Felkel).
Markers coloured by sigma (highest mean |r| with residual across expressions).

Per-expression performance:
  - Horton–Einstein: R²=0.87, RMSE=0.038
  - Pavlovskii: R²=0.85, RMSE=0.046
  - Lotter: R²=0.97, RMSE=0.011
  - Cox (arithmetic): R²=0.91, RMSE=0.030
  - Shear-force: R²=0.90, RMSE=0.033
  - Felkel: R²=0.95, RMSE=0.012

Residual structure is most associated with sigma.
### end fig5 ###

### fig6 ###
File     : fig6_pred_vs_obs_CF.png
Updated  : 2026-03-13 15:25
Notebook : notebooks/roughness_scale-cf-channel_eqs.ipynb
Saved in : ../figures/runaround_smooth
────────────────────────────────────────────────────────────────────────
Predicted vs observed n_e/<n> for four CF-field methods (Leading term, T0, Lotter, Felkel).
Markers coloured by sigma (highest |r| with leading-term residual).
RMSE/R² annotated per panel.

Per-method performance:
  - Leading term (CF): R²=0.56, RMSE=0.096
  - 2nd-order $T_0$ (CF): R²=0.13, RMSE=0.207
  - Lotter (CF): R² undefined, RMSE=0.239
  - Felkel: R²=0.51, RMSE=0.080

Residuals correlate most with sigma; dominant errors occur at extreme sigma values.
### end fig6 ###

### SI1 ###
File     : SI1_CF_errors_3panel.png
Updated  : 2026-03-13 15:25
Notebook : notebooks/roughness_scale-cf-channel_eqs.ipynb
Saved in : ../figures/runaround_smooth
────────────────────────────────────────────────────────────────────────
Three-panel correction-factor diagnostic.
Panels 1 & 2: relative U and h errors vs CF values, coloured by tr.
Panel 3: h/U identity comparison (correction factor only).
### end SI1 ###

### SI2 ###
File     : SI2_effect_ratio_vs_sigma_grid.png
Updated  : 2026-03-27 16:18
Notebook : notebooks/roughness_scale-analysis.ipynb
Saved in : ../figures/runaround_smooth
────────────────────────────────────────────────────────────────────────
$n_e/\\langle n \\rangle$ vs $\\sigma$ for all simulations at $f_V \\approx 0.6$, organised in a grid of panels by rainfall intensity ($p$, rows) and storm duration ($t_r$, columns). Points are coloured and shaped by vegetation arrangement type (or anisotropy when only blob patterns are present). Orientation categories are offset horizontally (±0.15 in $\\sigma$) for clarity.
### end SI2 ###
"""

# ---- Concise registry ----
concise_reg = f"""\
Figure Registry (concise)  •  {now}
Source: notebooks/roughness_scale-pattern.ipynb  |  Figures: /Users/octaviacrompton/Projects/roughness-scale/figures/runaround_smooth
========================================================================

Fig 1 — fig1_conceptual_sketch.png
  Conceptual sketch illustrating the equivalent-roughness problem on a heterogeneous hillslope ($f_V=0.4$, $n_{{\\rm veg}}=0.20$, $n_{{\\rm bare}}=0.05$). Panel (a) shows the patchy surface, (b) the uniform-roughness stack, and (c) the matched hydrograph, making the question '$n_e = ?$' explicit.

Fig 2 — fig2_effect_ratio_fV_2x2.png
  Main results figure: $n_e/\\langle n\\rangle$ vs $f_V$ for least and most intense storms (top row) with corresponding hydrographs at $f_V\\approx 0.5$ (bottom row), both coloured by patch lengthscale $L_V$. Demonstrates that the ratio is always below unity, is minimised near $f_V=0.5$, and increases with patch size.

Fig 3 — fig3_effect_ratio_vs_LV.png
  Effect ratio vs patch lengthscale $L_V$, panelled by $f_V$, coloured by orientation. Ratio increases with patch size toward unity; orientation sensitivity peaks near $f_V \\approx 0.5$.

Fig 4 — fig4_Sf_decomp_3x3.png
  $S_f$ decomposition figure (3×3): row 1 compares observed, $T_2$, and Lotter effect ratios vs $f_V$; rows 2–3 show the six dispersive-flux terms. Coloured by patch lengthscale $L_V$. Rows 2–3 share a y-axis so term magnitudes can be compared directly.

Fig 5 — fig5_obs_vs_pred_re_6panel.png
  Observed vs predicted r_e for six composite-roughness expressions, coloured by sigma. Lotter tracks the 1:1 line best (R²=0.97); Pavlovskii is weakest (R²=0.85).

Fig 6 — fig6_pred_vs_obs_CF.png
  Predicted vs observed n_e/<n> for four CF methods, coloured by sigma. Leading term (CF) performs best (R²=0.56); 2nd-order $T_0$ (CF) is weakest (R²=0.13).

SI1 — SI1_CF_errors_3panel.png
  CF diagnostic coloured by tr: relative U and h errors vs CF values, plus h/U identity for correction factor.

SI2 — SI2_effect_ratio_vs_sigma_grid.png
  Effect ratio vs sigma at fV≈0.6, panelled by p×tr, coloured by veg type / anisotropy. Reveals how spatial heterogeneity (sigma) and storm forcing jointly control the effective-roughness enhancement.
"""

(FIG_DIR / "figure_registry.txt").write_text(full_reg)
(FIG_DIR / "figure_registry_concise.txt").write_text(concise_reg)
print("  registries rewritten")

print("\nDone. Summary:")
print("  - fig3→fig4 (Sf_decomp, decomp notebook)")
print("  - fig4→fig5 (obs_vs_pred, cf notebook)")
print("  - fig5→fig6 (pred_vs_obs_CF, cf-channel_eqs notebook)")
print("  - NEW fig3 (effect_ratio_vs_LV, pattern notebook)")
print("  - PNGs renamed, registries rewritten")
