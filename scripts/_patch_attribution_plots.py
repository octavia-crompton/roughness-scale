#!/usr/bin/env python3
"""Patch cell 13 of decomp notebook: replace stacked bars with heatmap + lollipop,
soften palette, use _soft_colors for profile plot."""
import json, pathlib

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"
nb = json.loads(NB.read_text())

TARGET_MARKER = "Length-scale attribution of S_f decomposition terms"

for ci, cell in enumerate(nb["cells"]):
    if cell["cell_type"] != "code":
        continue
    src = "".join(cell["source"])
    if TARGET_MARKER not in src:
        continue

    # ---- 1. Replace _bar_colors + stacked bar figure with heatmap ----
    old_stacked = (
        "_x = np.arange(len(_plot_df))\n"
        "_bar_colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#d62728', '#9467bd', '#8c564b']\n"
        "\n"
        "fig_attr, ax_attr = plt.subplots(figsize=(10, 4))\n"
        "_bot_pos = np.zeros(len(_x))\n"
        "_bot_neg = np.zeros(len(_x))\n"
        "for _k, _lab in enumerate(_labels):\n"
        "    _vals = _plot_df[_lab].to_numpy(float)\n"
        "    _pos = np.where(_vals >= 0, _vals, 0)\n"
        "    _neg = np.where(_vals < 0, _vals, 0)\n"
        "    ax_attr.bar(_x, _pos, bottom=_bot_pos, color=_bar_colors[_k], label=_lab,\n"
        "                width=0.7, edgecolor='w', linewidth=0.3)\n"
        "    ax_attr.bar(_x, _neg, bottom=_bot_neg, color=_bar_colors[_k],\n"
        "                width=0.7, edgecolor='w', linewidth=0.3)\n"
        "    _bot_pos += _pos\n"
        "    _bot_neg += _neg\n"
        "\n"
        "_xlabels = [f'{v:g}' if isinstance(v, float) else v for v in _plot_df['fV']]\n"
        "ax_attr.set_xticks(_x)\n"
        "ax_attr.set_xticklabels(_xlabels, fontsize=FS_TICK)\n"
        "ax_attr.axhline(1, color='k', ls=':', lw=0.6)\n"
        "ax_attr.axhline(0, color='k', ls='-', lw=0.5)\n"
        "ax_attr.axvline(len(_x) - 1.5, color='grey', ls='--', lw=0.8)\n"
        "ax_attr.set_xlabel(r'$f_V$', fontsize=FS_LABEL)\n"
        "ax_attr.set_ylabel('fraction of Var$(S)$', fontsize=FS_LABEL)\n"
        "ax_attr.set_title('Covariance attribution of length-scale variance', fontsize=FS_TITLE)\n"
        "ax_attr.legend(fontsize=FS_LEG - 1, ncol=3, loc='upper left', frameon=True)\n"
        "ax_attr.tick_params(labelsize=FS_TICK)\n"
        "fig_attr.tight_layout()\n"
        "plt.show()"
    )

    new_heatmap = (
        "# -- Soft, muted colour palette for all figures in this cell --\n"
        "_soft_colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']\n"
        "\n"
        "# -- Figure: annotated heatmap of covariance attribution --\n"
        "_heat_data = _plot_df.set_index('fV')[_labels].astype(float)\n"
        "_nrows_h, _ncols_h = _heat_data.shape\n"
        "\n"
        "fig_attr, ax_attr = plt.subplots(figsize=(10, max(3.2, 0.55 * _nrows_h + 1.2)))\n"
        "_im = ax_attr.imshow(_heat_data.values, cmap='RdBu_r', vmin=-0.6, vmax=1.0,\n"
        "                     aspect='auto', interpolation='nearest')\n"
        "for _i in range(_nrows_h):\n"
        "    for _j in range(_ncols_h):\n"
        "        _v = _heat_data.iloc[_i, _j]\n"
        "        _tc = 'w' if abs(_v) > 0.35 else '0.15'\n"
        "        ax_attr.text(_j, _i, f'{_v:+.0%}', ha='center', va='center',\n"
        "                     fontsize=FS_TICK, fontweight='bold', color=_tc)\n"
        "# separator before avg row\n"
        "ax_attr.axhline(_nrows_h - 1.5, color='k', lw=1.2)\n"
        "ax_attr.set_xticks(range(_ncols_h))\n"
        "ax_attr.set_xticklabels(_labels, fontsize=FS_TICK - 1, rotation=25, ha='right')\n"
        "_ylabels_h = [f'$f_V={v:g}$' if isinstance(v, float) else v for v in _heat_data.index]\n"
        "ax_attr.set_yticks(range(_nrows_h))\n"
        "ax_attr.set_yticklabels(_ylabels_h, fontsize=FS_TICK)\n"
        "_cb = plt.colorbar(_im, ax=ax_attr, shrink=0.8, pad=0.02)\n"
        "_cb.set_label('fraction of Var$(S)$', fontsize=FS_LABEL)\n"
        "_cb.ax.tick_params(labelsize=FS_TICK)\n"
        "ax_attr.set_title('Covariance attribution of length-scale variance', fontsize=FS_TITLE)\n"
        "fig_attr.tight_layout()\n"
        "plt.show()"
    )

    if old_stacked not in src:
        print("ERROR: Could not find stacked bar block")
        exit(1)
    src = src.replace(old_stacked, new_heatmap)

    # ---- 2. Replace vertical bar slopes with horizontal lollipop ----
    old_slope = (
        "# -- Figure: grouped bar of slopes --\n"
        "fig_slope, ax_slope = plt.subplots(figsize=(8, 4))\n"
        "_sx = np.arange(_n_terms)\n"
        "_slope_vals = [_slopes[_lab] for _lab in _labels]\n"
        "_slope_ses  = [_slopes_se[_lab] for _lab in _labels]\n"
        "ax_slope.bar(_sx, _slope_vals, color=_bar_colors, width=0.6,\n"
        "             edgecolor='k', linewidth=0.4)\n"
        "ax_slope.errorbar(_sx, _slope_vals, yerr=_slope_ses, fmt='none',\n"
        "                  ecolor='k', capsize=3, lw=1)\n"
        "ax_slope.axhline(0, color='k', lw=0.5)\n"
        "ax_slope.axhline(_slope_S, color='grey', ls='--', lw=1,\n"
        "                 label=rf'total $\\partial S/\\partial L_V = {_slope_S:.4f}$')\n"
        "ax_slope.set_xticks(_sx)\n"
        "ax_slope.set_xticklabels(_labels, fontsize=FS_TICK - 1)\n"
        "ax_slope.set_ylabel(r'$\\partial T_k\\,/\\,\\partial L_V$  (m$^{-1}$)', fontsize=FS_LABEL)\n"
        "ax_slope.set_title(r'Sensitivity of each term to patch scale $L_V$ (fV fixed effects)',\n"
        "                   fontsize=FS_TITLE)\n"
        "ax_slope.legend(fontsize=FS_LEG, frameon=True)\n"
        "ax_slope.tick_params(labelsize=FS_TICK)\n"
        "fig_slope.tight_layout()\n"
        "plt.show()"
    )

    new_lollipop = (
        "# -- Figure: horizontal lollipop chart of slopes --\n"
        "fig_slope, ax_slope = plt.subplots(figsize=(8, 4))\n"
        "_sy = np.arange(_n_terms)\n"
        "_slope_vals = [_slopes[_lab] for _lab in _labels]\n"
        "_slope_ses  = [_slopes_se[_lab] for _lab in _labels]\n"
        "ax_slope.hlines(_sy, 0, _slope_vals, color=_soft_colors[:_n_terms], lw=2.5, zorder=2)\n"
        "ax_slope.scatter(_slope_vals, _sy, color=_soft_colors[:_n_terms], s=90, zorder=3,\n"
        "                 edgecolor='0.3', linewidth=0.5)\n"
        "ax_slope.errorbar(_slope_vals, _sy, xerr=_slope_ses, fmt='none',\n"
        "                  ecolor='0.4', capsize=3, lw=1, zorder=1)\n"
        "ax_slope.axvline(0, color='k', lw=0.5)\n"
        "ax_slope.axvline(_slope_S, color='0.5', ls='--', lw=1,\n"
        "                 label=rf'total $\\partial S/\\partial L_V = {_slope_S:.4f}$')\n"
        "ax_slope.set_yticks(_sy)\n"
        "ax_slope.set_yticklabels(_labels, fontsize=FS_TICK)\n"
        "ax_slope.set_xlabel(r'$\\partial T_k\\,/\\,\\partial L_V$  (m$^{-1}$)', fontsize=FS_LABEL)\n"
        "ax_slope.set_title(r'Sensitivity of each term to patch scale $L_V$ (fV fixed effects)',\n"
        "                   fontsize=FS_TITLE)\n"
        "ax_slope.legend(fontsize=FS_LEG, frameon=True, loc='lower right')\n"
        "ax_slope.tick_params(labelsize=FS_TICK)\n"
        "fig_slope.tight_layout()\n"
        "plt.show()"
    )

    if old_slope not in src:
        print("ERROR: Could not find slope bar block")
        exit(1)
    src = src.replace(old_slope, new_lollipop)

    # ---- 3. Soften profile plot colours ----
    src = src.replace("color=_bar_colors[_k]", "color=_soft_colors[_k]")

    # Write back
    cell["source"] = src.splitlines(keepends=True)
    print(f"Patched cell index {ci}")
    break
else:
    print("ERROR: target cell not found")
    exit(1)

NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
print("Done.")
