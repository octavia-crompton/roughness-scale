"""
1. cf-channel_eqs: append Lotter ≡ 1 proof to channel-equations markdown (cell 13)
2. cf-channel_eqs: promote SI2 → fig5 in save cell (cell 19)
3. notebook 4:     remove channel-equation cells 13-18
"""
import json
from pathlib import Path

NB_CH  = Path("notebooks/roughness_scale-cf-channel_eqs.ipynb")
NB_CF  = Path("notebooks/4. roughness_scale-cf.ipynb")

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  cf-channel_eqs: add Lotter ≡ 1 proof to channel-equations markdown
# ═══════════════════════════════════════════════════════════════════════════════
with open(NB_CH) as f:
    nb_ch = json.load(f)

LOTTER_PROOF = r"""
---

### Why Lotter $\equiv 1$ for CF inputs

The CF depth is $h_i = (q_0 l n_i / S_0^{1/2})^{3/5}$, so raising to the $5/3$ power gives

$$h_i^{5/3} = \frac{q_0\, l}{S_0^{1/2}}\, n_i \equiv C\, n_i$$

where $C = q_0 l / S_0^{1/2}$ is the same for every cell. Therefore $\langle h^{5/3}\rangle = C\langle n\rangle$ and $h_i^{5/3}/n_i = C$ (a cell-independent constant), so $\langle h^{5/3}/n\rangle = C$. Substituting:

$$\frac{n_e}{\langle n\rangle} = \frac{C\langle n\rangle}{\langle n\rangle \cdot C} = 1$$

The Lotter formula is uninformative for any set of inputs in which $h^{5/3} \propto n$ — which is exactly what the CF approximation enforces.
"""

cell_13 = nb_ch["cells"][13]
src_13 = "".join(cell_13.get("source", []))
if "Lotter" in src_13 and "uninformative" in src_13:
    print("  Lotter proof already present in cell 13 — skipping")
else:
    # Append proof to end of existing markdown
    old_lines = cell_13["source"]
    # Make sure last line has a newline
    if old_lines and not old_lines[-1].endswith("\n"):
        old_lines[-1] += "\n"
    cell_13["source"] = old_lines + LOTTER_PROOF.splitlines(keepends=True)
    print("  Added Lotter ≡ 1 proof to cell 13")

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  cf-channel_eqs: promote SI2 → fig5 in save cell (cell 19)
# ═══════════════════════════════════════════════════════════════════════════════
cell_19 = nb_ch["cells"][19]
src_19_lines = cell_19["source"]
new_lines_19 = []
for line in src_19_lines:
    line = line.replace("_si2_name", "_fig5_name")
    line = line.replace("SI2_channel_eqs_9panel.png", "fig5_obs_vs_pred_re_9panel.png")
    line = line.replace("'SI2'", "'fig5'")
    line = line.replace("Save absolute comparison as SI2",
                        "Save absolute comparison as fig5")
    line = line.replace("Saved  {_si2_name}", "Saved  {_fig5_name}")
    new_lines_19.append(line)
cell_19["source"] = new_lines_19
cell_19["outputs"] = []
print("  Promoted SI2 → fig5 in cell 19")

with open(NB_CH, "w") as f:
    json.dump(nb_ch, f, indent=1, ensure_ascii=False)
    f.write("\n")
print("  Wrote cf-channel_eqs.ipynb")

# ═══════════════════════════════════════════════════════════════════════════════
# 3.  notebook 4: remove channel-equation cells 13-18
# ═══════════════════════════════════════════════════════════════════════════════
with open(NB_CF) as f:
    nb_cf = json.load(f)

# Verify expected cells before deleting
cell13_src = "".join(nb_cf["cells"][13].get("source", []))
cell18_src = "".join(nb_cf["cells"][18].get("source", []))
assert "Channel equation" in cell13_src, f"Cell 13 unexpected: {cell13_src[:60]}"
assert "fig5" in cell18_src or "n_exp" in cell18_src, f"Cell 18 unexpected: {cell18_src[:60]}"

removed = nb_cf["cells"][13:19]
for i, c in enumerate(removed):
    src = "".join(c.get("source", []))[:60].replace("\n", " | ")
    print(f"  Removing cell [{13+i}]: {src}")

nb_cf["cells"] = nb_cf["cells"][:13] + nb_cf["cells"][19:]
print(f"  Notebook 4 now has {len(nb_cf['cells'])} cells (was {len(nb_cf['cells'])+6})")

with open(NB_CF, "w") as f:
    json.dump(nb_cf, f, indent=1, ensure_ascii=False)
    f.write("\n")
print("  Wrote 4. roughness_scale-cf.ipynb")
print("\nDone.")
