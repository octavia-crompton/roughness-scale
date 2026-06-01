"""
Fix the cf-channel_eqs notebook:
  1. Rewrite markdown cell to clearly separate direct CF vs perturbation expansion
  2. Replace T₁ with full 2nd-order in the plot cell
"""
import json
from pathlib import Path

NB_PATH = Path("notebooks/roughness_scale-cf-channel_eqs.ipynb")

with open(NB_PATH) as f:
    nb = json.load(f)

# ═══════════════════════════════════════════════════════════════════════════════
# 1.  Rewrite the markdown cell
# ═══════════════════════════════════════════════════════════════════════════════
NEW_MARKDOWN = r"""### Two approaches to predicting $n_e/\langle n\rangle$ from the roughness field

Both methods use only the roughness field $n(x,y)$ — no simulation output — with the kinematic-wave CF approximation providing per-cell depth $h_i \propto n_i^{3/5}$ and velocity $U_i \propto n_i^{-3/5}$.

---

#### Approach 1 — Direct back-calculation (Leading term)

Spatially average the CF-estimated depth and velocity, then back-calculate $n_e$ from Manning's equation:

$$\frac{n_e}{\langle n\rangle} = \frac{\langle h\rangle_{\rm CF}^{2/3}\, S_0^{1/2}}{\langle U\rangle_{\rm CF}\,\langle n\rangle}$$

This is **exact** given the CF inputs — no expansion or truncation. It preserves the full nonlinearity of the spatial average.

---

#### Approach 2 — Perturbation expansion of $\langle S_f \rangle$

Write each field as a mean plus normalised fluctuation: $n = \bar{n}(1+r)$, $U = \bar{U}(1+\upsilon)$, $h = \bar{h}(1+\eta)$. Setting $\langle S_f \rangle = n_e^2 \bar{U}^2 \bar{h}^{-4/3}$ gives

$$(n_e/\bar{n})^2 = \bigl\langle (1+r)^2 (1+\upsilon)^2 (1+\eta)^{-4/3} \bigr\rangle$$

Expanding to **second order** (first-order means vanish):

$$\frac{n_e}{\langle n\rangle} \approx 1 + \underbrace{\langle r^2\rangle + \langle \upsilon^2\rangle + 4\langle r\upsilon\rangle}_{T_0:\; r,\,\upsilon\text{ terms only}} \underbrace{- \tfrac{8}{3}\langle r\eta\rangle - \tfrac{8}{3}\langle \eta\upsilon\rangle + \tfrac{14}{9}\langle \eta^2\rangle}_{\text{depth-coupling terms}}$$

All six covariance terms are second-order. The $T_0$ / full-2nd-order split separates terms that involve only roughness and velocity from those that also involve depth.

**After CF substitution** ($\upsilon \approx -\tfrac{3}{5}r$, $\eta \approx \tfrac{3}{5}r$), every term becomes a multiple of $\langle r^2\rangle$:

| Term | CF substitution | Coefficient on $\langle r^2\rangle$ |
|------|-----------------|--------------------------------------|
| $\langle r^2\rangle$ | — | $+1$ |
| $\langle \upsilon^2\rangle$ | $\tfrac{9}{25}\langle r^2\rangle$ | $+\tfrac{9}{25}$ |
| $4\langle r\upsilon\rangle$ | $-\tfrac{12}{5}\langle r^2\rangle$ | $-\tfrac{12}{5}$ |
| $-\tfrac{8}{3}\langle r\eta\rangle$ | $-\tfrac{8}{5}\langle r^2\rangle$ | $-\tfrac{8}{5}$ |
| $-\tfrac{8}{3}\langle \eta\upsilon\rangle$ | $+\tfrac{24}{25}\langle r^2\rangle$ | $+\tfrac{24}{25}$ |
| $\tfrac{14}{9}\langle \eta^2\rangle$ | $\tfrac{14}{25}\langle r^2\rangle$ | $+\tfrac{14}{25}$ |

Summing:

| Level | CF formula |
|-------|------------|
| **$T_0$ (CF)** — $r, \upsilon$ only | $n_e/\langle n\rangle \approx 1 - \tfrac{26}{25}\langle r^2\rangle$ |
| **Full 2nd-order (CF)** — all six terms | $n_e/\langle n\rangle \approx 1 - \tfrac{28}{25}\langle r^2\rangle$ |

The depth-coupling terms nearly cancel under CF ($-\frac{8}{5} + \frac{24}{25} + \frac{14}{25} = -\frac{2}{25}$), adding only a small correction beyond $T_0$. The perturbation expansion should underperform the direct back-calculation because truncation at second order loses higher-order contributions.
"""

for i, cell in enumerate(nb["cells"]):
    if cell.get("cell_type") != "markdown":
        continue
    src = "".join(cell.get("source", []))
    if "Perturbation expansion" in src and "CF substitution" in src:
        cell["source"] = NEW_MARKDOWN.splitlines(keepends=True)
        cell["outputs"] = []
        print(f"  Replaced markdown cell at index {i}")
        break
else:
    print("  WARNING: markdown cell not found")

# ═══════════════════════════════════════════════════════════════════════════════
# 2.  In the plot cell: replace er_T1_CF with er_T2_CF (full 2nd-order)
# ═══════════════════════════════════════════════════════════════════════════════
for i, cell in enumerate(nb["cells"]):
    src_lines = cell.get("source", [])
    src = "".join(src_lines)
    if "er_T1_CF" in src and "methods = [" in src:
        new_lines = []
        for line in src_lines:
            # Replace T1 computation with full T2 (all six terms)
            if line.strip().startswith("# T") and "adds" in line and "cross-terms" in line:
                new_lines.append(
                    "# Full 2nd-order: all six terms  (CF substitutions)\n"
                )
            elif "er_T1_CF" in line and "summary[" in line and "=" in line and "(" in line:
                new_lines.append(
                    "summary['er_T2_CF'] = (\n"
                )
            elif line.strip().startswith("- (8/3) * (3/5)"):
                new_lines.append(line)  # keep -8/3 <r eta>
            elif line.strip().startswith("- (8/3) * (-9/25)"):
                new_lines.append(line)  # keep -8/3 <eta ups>
                # Add T2 term after
                new_lines.append(
                    "    + (14/9) * (9/25) * summary['<r2>']       # 14/9 ⟨η²⟩\n"
                )
            elif "'er_T1_CF'" in line:
                new_lines.append(line.replace("er_T1_CF", "er_T2_CF").replace(
                    r"$T_1$: adds $\eta$ cross-terms (CF)",
                    r"Full 2nd-order (CF)"
                ))
            elif "er_T1_CF" in line:
                new_lines.append(line.replace("er_T1_CF", "er_T2_CF"))
            elif "three CF perturbation expansions" in line:
                new_lines.append(line.replace(
                    "three CF perturbation expansions",
                    "three CF-based predictions"
                ))
            elif "T0, T1" in line:
                new_lines.append(line.replace("T0, T1", "T0, full 2nd-order"))
            elif "T0 with r" in line and "T1 adding" in line:
                new_lines.append(line.replace(
                    "T0 with r/υ terms, T1 adding η cross-terms",
                    "T0 r/υ terms only, full 2nd-order all six terms"
                ))
            else:
                new_lines.append(line)
        cell["source"] = new_lines
        cell["outputs"] = []
        print(f"  Updated plot cell at index {i}")
        break
else:
    print("  WARNING: plot cell not found")

# ═══════════════════════════════════════════════════════════════════════════════
with open(NB_PATH, "w") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write("\n")
print("  Done")
