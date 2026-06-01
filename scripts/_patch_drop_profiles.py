#!/usr/bin/env python3
"""Remove the conditional profiles plot (plot 3) from the attribution cell."""
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

# Remove the entire conditional profiles section
PROFILES_BLOCK = """\

# ============================================================================
# 3. Conditional profiles -- E[T_k | sigma] across patch scales
# ============================================================================
_sigma_levels = sorted(_sub['sigma'].dropna().unique())
_lv_vals = [_sigma_to_LV.get(s, s) for s in _sigma_levels]

_profiles = {}
for _k, _lab in enumerate(_labels):
    _means = []
    for _s in _sigma_levels:
        _means.append(_sub.loc[_sub['sigma'] == _s, _Tcols[_k]].mean())
    _profiles[_lab] = _means

fig_prof, ax_prof = plt.subplots(figsize=(7, 4.5))
for _k, _lab in enumerate(_labels):
    ax_prof.plot(_lv_vals, _profiles[_lab], 'o-', color=_earth_colors[_k],
                 label=_lab, markersize=5, lw=1.5)
ax_prof.axhline(0, color='k', ls='-', lw=0.5)
ax_prof.set_xlabel(r'$L_V$ (m)', fontsize=FS_LABEL)
ax_prof.set_ylabel('mean term value', fontsize=FS_LABEL)
ax_prof.set_title(r'Conditional mean of each $S_f$ term vs patch scale',
                  fontsize=FS_TITLE)
ax_prof.legend(fontsize=FS_LEG - 1, frameon=True)
ax_prof.tick_params(labelsize=FS_TICK)
fig_prof.tight_layout()
plt.show()"""

src = src.replace(PROFILES_BLOCK, "")

# Update save section: 3 figures → 2 figures
src = src.replace(
    "# -- Save all three to scratch --\n"
    "_, _scratch, _ = _fig_dirs()\n"
    "fig_attr.savefig(_os.path.join(_scratch, 'lengthscale_cov_attribution.png'),\n"
    "                 dpi=200, bbox_inches='tight')\n"
    "fig_slope.savefig(_os.path.join(_scratch, 'lengthscale_slope_decomp.png'),\n"
    "                  dpi=200, bbox_inches='tight')\n"
    "fig_prof.savefig(_os.path.join(_scratch, 'lengthscale_conditional_profiles.png'),\n"
    "                 dpi=200, bbox_inches='tight')\n"
    "print(\"\\nSaved 3 figures to scratch/\")",

    "# -- Save to scratch --\n"
    "_, _scratch, _ = _fig_dirs()\n"
    "fig_attr.savefig(_os.path.join(_scratch, 'lengthscale_cov_attribution.png'),\n"
    "                 dpi=200, bbox_inches='tight')\n"
    "fig_slope.savefig(_os.path.join(_scratch, 'lengthscale_slope_decomp.png'),\n"
    "                  dpi=200, bbox_inches='tight')\n"
    "print(\"\\nSaved 2 figures to scratch/\")",
)

cell["source"] = src.splitlines(keepends=True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print("OK – removed conditional profiles plot")
