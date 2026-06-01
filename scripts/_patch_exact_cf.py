#!/usr/bin/env python3
"""Patch notebook 4: replace linearized CF substitutions with exact (nonlinear)
CF moments in both the markdown explanation and the code cell."""
import json, pathlib, re

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "4. roughness_scale-cf-channel_eqs.ipynb"
nb = json.loads(NB.read_text())
cells = nb["cells"]

changes = 0

# ═══════════════════════════════════════════════════════════════════════════════
# 1. MARKDOWN CELL — replace the "After CF substitution" section
# ═══════════════════════════════════════════════════════════════════════════════
NEW_MD_SECTION = r"""**Exact CF moments** — Under the CF power law, $U_i \propto n_i^{-3/5}$ and $h_i \propto n_i^{3/5}$, so the normalised fluctuations are

$$1+\upsilon_i = (1+r_i)^{-3/5}, \qquad 1+\eta_i = (1+r_i)^{3/5}$$

exactly (no linearisation). For a binary roughness field ($n_v$, $n_b$ with fraction $f_V$) we compute six second moments directly from these:

$$\langle r^2\rangle,\; \langle\upsilon^2_{\rm CF}\rangle,\; \langle r\,\upsilon_{\rm CF}\rangle,\; \langle r\,\eta_{\rm CF}\rangle,\; \langle\eta_{\rm CF}\upsilon_{\rm CF}\rangle,\; \langle\eta^2_{\rm CF}\rangle$$

and plug them into the second-order expansion. Because the exact CF fluctuations include all orders of $r$, the resulting moments are more accurate than the linearised substitutions ($\upsilon \approx -\tfrac{3}{5}r$, $\eta \approx \tfrac{3}{5}r$), which collapse every term to a multiple of $\langle r^2\rangle$ and introduce systematic bias.

> **Note:** Under exact CF, the pointwise product $(1+r)^2(1+\upsilon)^2(1+\eta)^{-4/3} = S_0/(\bar{n}^2\bar{U}^2\bar{h}^{-4/3})$ is the same constant for every cell, so the perturbation formula with exact CF moments recovers the leading-term result identically. The T₀ / full split then isolates how much of the Jensen's-inequality effect is captured by the $r$–$\upsilon$ terms alone versus the depth-coupling terms.
"""

for ci, cell in enumerate(cells):
    if cell["cell_type"] != "markdown":
        continue
    src = "".join(cell["source"])
    if "After CF substitution" in src and "Summing:" in src:
        # Find everything from "**After CF substitution**" to end of cell
        lines = cell["source"]
        new_lines = []
        found_start = False
        for line in lines:
            if "**After CF substitution**" in line:
                found_start = True
                # Insert new section instead
                new_lines.append(NEW_MD_SECTION)
                break
            new_lines.append(line)
        if found_start:
            cell["source"] = new_lines
            changes += 1
            print(f"  Updated markdown cell {ci}: replaced linearised CF table with exact CF explanation")
        break
else:
    print("WARNING: could not find markdown cell with 'After CF substitution'")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. CODE CELL — replace linearised T₀/T₂ with exact CF moments
# ═══════════════════════════════════════════════════════════════════════════════
NEW_CODE_BLOCK = [
    "# T₀ and T₂ using EXACT CF moments (no linearisation)\n",
    "# Under CF: U_i ∝ n_i^{-3/5}, h_i ∝ n_i^{3/5}\n",
    "# So (1+υ) = (1+r)^{-3/5} and (1+η) = (1+r)^{3/5} exactly.\n",
    "_nbar  = summary['<n>']\n",
    "_rv    = (summary['alpha_v'] - _nbar) / _nbar   # r for veg cells\n",
    "_rb    = (summary['alpha_b'] - _nbar) / _nbar   # r for bare cells\n",
    "_fV    = summary['fV']\n",
    "\n",
    "# exact CF fluctuations (not linearised)\n",
    "_uv_v  = (1 + _rv)**(-3/5) - 1    # υ for veg cells\n",
    "_uv_b  = (1 + _rb)**(-3/5) - 1    # υ for bare cells\n",
    "_eta_v = (1 + _rv)**(3/5)  - 1    # η for veg cells\n",
    "_eta_b = (1 + _rb)**(3/5)  - 1    # η for bare cells\n",
    "\n",
    "# six second moments (binary-weighted)\n",
    "_r2_cf      = _fV * _rv**2       + (1-_fV) * _rb**2\n",
    "_ups2_cf    = _fV * _uv_v**2     + (1-_fV) * _uv_b**2\n",
    "_r_ups_cf   = _fV * _rv*_uv_v    + (1-_fV) * _rb*_uv_b\n",
    "_r_eta_cf   = _fV * _rv*_eta_v   + (1-_fV) * _rb*_eta_b\n",
    "_eta_ups_cf = _fV * _eta_v*_uv_v + (1-_fV) * _eta_b*_uv_b\n",
    "_eta2_cf    = _fV * _eta_v**2    + (1-_fV) * _eta_b**2\n",
    "\n",
    "# T₀: r, υ terms only\n",
    "summary['er_T0_CF'] = 1 + _r2_cf + _ups2_cf + 4*_r_ups_cf\n",
    "\n",
    "# Full 2nd-order: all six terms\n",
    "summary['er_T2_CF'] = (\n",
    "    summary['er_T0_CF']\n",
    "    - (8/3) * _r_eta_cf\n",
    "    - (8/3) * _eta_ups_cf\n",
    "    + (14/9) * _eta2_cf\n",
    ")\n",
]

for ci, cell in enumerate(cells):
    if cell["cell_type"] != "code":
        continue
    src = "".join(cell["source"])
    if "er_T0_CF" in src and "er_T2_CF" in src and "er_lead_CF" in src:
        lines = cell["source"]
        new_lines = []
        skip = False
        for line in lines:
            # Start skipping at the old T₀ computation
            if "T₀: keeps r, υ terms only" in line or "T₀:" in line and "linearised" in "".join(lines):
                skip = True
            if skip and "# T₀: keeps r" in line:
                skip = True
            if not skip:
                # Check for the old linearised block start
                if "# T₀: keeps r, υ terms only" in line:
                    skip = True
                    continue
                new_lines.append(line)
            else:
                # Stop skipping after the er_T2_CF block ends (blank line or next section)
                if "# ── define 3-panel" in line:
                    skip = False
                    new_lines.extend(NEW_CODE_BLOCK)
                    new_lines.append("\n")
                    new_lines.append(line)
        cell["source"] = new_lines
        changes += 1
        print(f"  Updated code cell {ci}: replaced linearised with exact CF moments")
        break
else:
    print("WARNING: could not find code cell with er_T0_CF / er_T2_CF")


if changes:
    NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
    print(f"\nDone — {changes} change(s) written.")
else:
    print("No changes needed.")
