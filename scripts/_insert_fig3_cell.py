"""Insert fig3 save cell into pattern notebook at position 22 (after LV figure cell)."""
import json

nb_path = "notebooks/2. roughness_scale-pattern.ipynb"
with open(nb_path, "r") as f:
    nb = json.load(f)

# Verify target cell
target_idx = 21  # JSON 0-indexed
assert nb["cells"][target_idx]["id"] == "9c4badb8", (
    f"unexpected cell id: {nb['cells'][target_idx].get('id')}")
print(f"Target cell: {nb['cells'][target_idx]['source'][0][:60]}...")

save_cell = {
    "cell_type": "code",
    "metadata": {},
    "outputs": [],
    "source": [
        "# \u2500\u2500 Save fig3 \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n",
        "_fig_dir, _, _ = _fig_dirs()\n",
        "_name = 'fig3_effect_ratio_vs_LV.png'\n",
        "fig.savefig(_os.path.join(_fig_dir, _name), dpi=300, bbox_inches='tight')\n",
        "update_figure_registry(\n",
        "    'fig3', _name,\n",
        "    description=(\n",
        "        '2\u00d73 panel figure showing $n_e/\\\\langle n\\\\rangle$ vs patch lengthscale $L_V$ '\n",
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

nb["cells"].insert(22, save_cell)
with open(nb_path, "w") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write("\n")
print(f"Inserted fig3 save cell at index 22 (after LV figure cell)")
