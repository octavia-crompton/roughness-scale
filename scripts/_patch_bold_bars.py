#!/usr/bin/env python3
"""Patch cell 13: bold palette + horizontal bar chart instead of lollipop."""
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

# ── 1. Replace pastel palette with bold, saturated colours ──
src = src.replace(
    "# -- Soft, muted colour palette for all figures in this cell --\n"
    "_soft_colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']",
    "# -- Bold, saturated colour palette for all figures in this cell --\n"
    "_bold_colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']",
)

# ── 2. Replace lollipop chart with horizontal bar chart ──
OLD_LOLLIPOP = """\
# -- Figure: horizontal lollipop chart of slopes --
fig_slope, ax_slope = plt.subplots(figsize=(8, 4))
_sy = np.arange(_n_terms)
_slope_vals = [_slopes[_lab] for _lab in _labels]
_slope_ses  = [_slopes_se[_lab] for _lab in _labels]
ax_slope.hlines(_sy, 0, _slope_vals, color=_soft_colors[:_n_terms], lw=2.5, zorder=2)
ax_slope.scatter(_slope_vals, _sy, color=_soft_colors[:_n_terms], s=90, zorder=3,
                 edgecolor='0.3', linewidth=0.5)
ax_slope.errorbar(_slope_vals, _sy, xerr=_slope_ses, fmt='none',
                  ecolor='0.4', capsize=3, lw=1, zorder=1)
ax_slope.axvline(0, color='k', lw=0.5)
ax_slope.axvline(_slope_S, color='0.5', ls='--', lw=1,
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

NEW_BAR = """\
# -- Figure: horizontal bar chart of slopes --
fig_slope, ax_slope = plt.subplots(figsize=(8, 4))
_sy = np.arange(_n_terms)
_slope_vals = [_slopes[_lab] for _lab in _labels]
_slope_ses  = [_slopes_se[_lab] for _lab in _labels]
ax_slope.barh(_sy, _slope_vals, xerr=_slope_ses, height=0.6,
              color=_bold_colors[:_n_terms], edgecolor='0.2', linewidth=0.6,
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

src = src.replace(OLD_LOLLIPOP, NEW_BAR)

# ── 3. Update all remaining _soft_colors refs to _bold_colors ──
src = src.replace("_soft_colors", "_bold_colors")

# Write back
cell["source"] = src.splitlines(keepends=True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print("OK – patched palette + bar chart")
