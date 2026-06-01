#!/usr/bin/env python3
"""Remove redundancies in overleaf/main.tex as identified in the audit.

Changes:
1. Delete the one-line "Channel flow" stub in the Introduction (lines 71-73)
2. Trim the redundant research-question framing at end of Introduction
   (merge into Problem statement)
3. In the appendix "Discharge estimate with h": remove the re-derivation of
   n_e/⟨n⟩ = ⟨h⟩^{5/3} / (⟨n⟩ ⟨h^{5/3}/n⟩) and "This is the same as the
   Lotter equation" — keep only the normalized-fluctuation expansion (the new part).
4. Remove duplicate "Expand 1/n" and "Expand (h/⟨h⟩)^{5/3}" paragraphs
   (already stated in the "Second-order approximation" preamble just above).
5. Merge "Definitions" + "Force Balance" + "Connecting Lotter equation to shear
   stress" into one clean appendix subsection.
"""
import pathlib, sys

tex_path = pathlib.Path(__file__).resolve().parent.parent / "overleaf" / "main.tex"
tex = tex_path.read_text()
original = tex  # keep for comparison

# ═══════════════════════════════════════════════════════════════════════════
# 1. Delete the Channel flow stub (lines 71-73)
# ═══════════════════════════════════════════════════════════════════════════
old_channel = (
    "\n\\subsubsection*{Channel flow}\n"
    "\n"
    "In 1993, Lotter derived the $n_e$ value based on the assumption that "
    "total discharge equals the sum of discharges of each sub-section.\n"
)
assert old_channel in tex, "Channel flow stub not found"
tex = tex.replace(old_channel, "\n")

# ═══════════════════════════════════════════════════════════════════════════
# 2. Trim redundant research question at end of Introduction.
#    The full "central question" sentence + the Problem statement opening
#    overlap. Remove the last sentence of the Introduction paragraph and
#    the first paragraph of Problem statement, replacing with a cleaner
#    bridge.
# ═══════════════════════════════════════════════════════════════════════════
old_bridge = (
    "Despite this body of work, the problem of how spatial pattern influences "
    "calibration roughness remains understudied. Relatively few studies have "
    "explicitly connected the spatial pattern of vegetation roughness — cover "
    "fraction, patch size, anisotropy, connectivity — to the equivalent roughness "
    "that a calibrated uniform model would recover.  This is the central question "
    "of the present study: we relate spatial heterogeneity in surface roughness "
    "(e.g., vegetation patches) to an \\emph{equivalent} (i.e., calibration, "
    "homogeneous) roughness parameter that reproduces a chosen bulk response "
    "(e.g., outlet hydrograph, domain-mean depth, infiltration fraction).\n"
    "\n"
    "\n"
    "\\section{Problem statement}\n"
    "\n"
    "Overland flow resistance is challenging to predict because it depends on the "
    "spatial pattern of roughness elements, which can vary widely across landscapes "
    "and influence flow in complex ways. \n"
    "The equivalent roughness parameter that emerges from calibration to bulk "
    "responses (e.g., runoff coefficient, infiltration fraction) can differ "
    "substantially from the spatial average of local roughness values, and this "
    "discrepancy depends on the spatial arrangement of roughness elements. "
    "Approaches to average spatially distributed roughness include empirical "
    "regression models \\citep{kim_hydraulic_2012}, analytical upscaling of the "
    "Navier--Stokes equations \\citep{nikora_double-averaging_2007}, and "
    "numerical simulations of flow over heterogeneous surfaces "
    "\\citep{cea_experimental_2014}, but a systematic understanding of how "
    "spatial pattern — patch size, connectivity, and anisotropy — influences "
    "$n_e$ remains lacking.\n"
    "\n"
    "In composite channels, $n_e$ falls below the area-weighted mean because "
    "smooth sub-reaches carry disproportionately more discharge "
    "\\citep{djajadi2009comparative, flintham1992manning}; we ask whether an "
    "analogous mechanism operates on vegetated hillslopes, and what spatial "
    "pattern statistics control the magnitude of $n_e/\\langle n \\rangle$."
)
new_bridge = (
    "Despite this body of work, relatively few studies have explicitly connected "
    "the spatial pattern of vegetation roughness — cover fraction, patch size, "
    "anisotropy, connectivity — to the equivalent roughness that a calibrated "
    "uniform model would recover.\n"
    "\n"
    "\n"
    "\\section{Problem statement}\n"
    "\n"
    "The equivalent roughness $n_e$ that emerges from calibrating a uniform model "
    "to bulk responses (e.g., outlet hydrograph, infiltration fraction) can differ "
    "substantially from the spatial average $\\langle n \\rangle$. "
    "In composite channels, $n_e$ falls below the area-weighted mean because "
    "smooth sub-reaches carry disproportionately more discharge "
    "\\citep{djajadi2009comparative, flintham1992manning}; we ask whether an "
    "analogous mechanism operates on vegetated hillslopes, and what spatial "
    "pattern statistics control the magnitude of $n_e/\\langle n \\rangle$."
)
assert old_bridge in tex, "Research-question bridge not found"
tex = tex.replace(old_bridge, new_bridge)

# ═══════════════════════════════════════════════════════════════════════════
# 3. Consolidate "Discharge estimate with h": remove the re-derivation
#    of the Lotter equation (already in Discussion) and keep only the
#    normalized-fluctuation expansion.
# ═══════════════════════════════════════════════════════════════════════════
old_h_section = (
    "\\subsubsection*{Discharge estimate with $h$}\n"
    "\n"
    " Derivation and interpretation of why the covariance term "
    "$\\langle n' U' \\rangle$ predicts the effect ratio when $n_e$ "
    "is defined via hydrograph calibration.\n"
    "\n"
    "For depth\u2013averaged flow with Manning's equation,\n"
    "\n"
    "\\begin{align}\n"
    "U &= \\frac{1}{n}\\, h^{2/3} S_0^{1/2}, \\\\\n"
    "q &= U h \\;=\\; \\frac{1}{n}\\, h^{5/3} S_0^{1/2}.\n"
    "\\end{align}\n"
    "At a given time, the patchy hillslope carries\n"
    "\\begin{equation}\n"
    "Q_{\\text{het}} \\;=\\; \\iint q \\, dA \\;=\\; S_0^{1/2} "
    "\\iint \\frac{h^{5/3}}{n}\\, dA.\n"
    "\\end{equation}\n"
    "\n"
    "Define the equivalent uniform case (same mean depth "
    "$\\langle h\\rangle$ over area $A$) by\n"
    "\\begin{equation}\n"
    "Q_{\\text{uni}}(n_e) \\;=\\; A \\,\\frac{1}{n_e}\\, "
    "\\langle h\\rangle^{5/3} S_0^{1/2}.\n"
    "\\end{equation}\n"
    "\n"
    "Setting  $Q_{\\text{uni}}=Q_{\\text{het}}$ and simplifying:\n"
    "\n"
    "\\begin{equation}\n"
    "\\frac{1}{n_e}\\,\\langle h\\rangle^{5/3} \\;=\\; "
    "\\left\\langle \\frac{h^{5/3}}{n} \\right\\rangle\n"
    "\\end{equation}\n"
    "\\begin{equation}\n"
    "\\frac{n_e}{\\langle n\\rangle}\n"
    "=\n"
    "\\frac{\\langle h\\rangle^{5/3}}{\\langle n\\rangle \\,"
    "\\Big\\langle \\frac{h^{5/3}}{n} \\Big\\rangle }.\n"
    "\\end{equation}\n"
    "\n"
    "The equivalent roughness is less than the spatial average because\n"
    "$\\langle h\\rangle^{5/3} < \\langle h^{5/3} \\rangle$ . \n"
    "The ratio $\\frac{ \\langle h^{5/3}\\rangle }"
    "{  \\langle h\\rangle^{5/3} }$ increases with increasing "
    "vegetation patch scale, so equivalent roughness decreases.\n"
    "\n"
    "\n"
    "This is the same as the Lotter equation: \n"
    "\\begin{equation}\n"
    "{n_e}\n"
    "= \\langle h^{5/3} \\rangle \\bigg \\langle  \\frac{ h^{5/3}}"
    "{n}\\bigg \\rangle^{-1} \n"
    "\\end{equation}\n"
    "\n"
    "\n"
    "Define: \n"
)
new_h_section = (
    "\\subsubsection*{Second-order expansion in terms of $h$ and $n$}\n"
    "\n"
    "Starting from the equal-discharge result derived in the Discussion "
    "(Eq.~above), we expand in normalized fluctuations.\n"
    "\n"
    "Define: \n"
)
assert old_h_section in tex, "Discharge estimate with h section not found"
tex = tex.replace(old_h_section, new_h_section)

# ═══════════════════════════════════════════════════════════════════════════
# 4. Remove duplicate "Expand 1/n" and "Expand (h/⟨h⟩)^{5/3}" paragraphs.
#    The second-order preamble (lines 532-544) already states these
#    expansions concisely; the longer re-derivation that follows is
#    redundant.
# ═══════════════════════════════════════════════════════════════════════════
old_dup_expand = (
    "\\paragraph{Expand  $1/n$:} \n"
    "\n"
    "\\begin{align}\n"
    "\\frac{1}{n}\n"
    "= \\frac{1}{\\langle n\\rangle + n'}\n"
    "% = \\frac{1}{\\langle n\\rangle}\\,\\frac{1}{1 + \\tfrac{n'}"
    "{\\langle n\\rangle}}\n"
    "= \\frac{1}{\\langle n\\rangle}\\,(1+r)^{-1},\n"
    "\\quad\\text{where } r \\equiv \\frac{n'}{\\langle n\\rangle}.\n"
    "\\end{align}\n"
    "\n"
    "For $|r|<1$, the Taylor series is:\n"
    "\n"
    "\\begin{align}\n"
    "(1+r)^{-1} = \\sum_{k=0}^{\\infty}(-r)^k\n"
    "= 1 - r + r^2 - r^3 + \\cdots .\n"
    "\\end{align}\n"
    "\\begin{align}\n"
    "\\frac{1}{n}\n"
    "\\approx \\frac{1}{\\langle n\\rangle}\\!\\left(1 - "
    "\\frac{n'}{\\langle n\\rangle} + \\frac{n'^2}"
    "{\\langle n\\rangle^2}\\right)\n"
    "\\end{align}\n"
    "\n"
    "\\paragraph{Expand $\\big(h/\\langle h\\rangle\\big)^{5/3}$.}\n"
    "Write $h=\\langle h\\rangle(1+\\eta)$ with "
    "$\\eta=h'/\\langle h\\rangle$.\n"
    "\n"
    "\\begin{align}\n"
    "\\left(\\frac{h}{\\langle h\\rangle}\\right)^{5/3}\n"
    "= (1+\\eta)^{5/3}\n"
    "\\end{align}\n"
    "\n"
    "\n"
    "\\begin{align}\n"
    "\\left(\\frac{h}{\\langle h\\rangle}\\right)^{5/3}\n"
    "\\approx 1 + \\frac{5}{3}\\eta + \\frac{5}{9}\\eta^2\n"
    "\\end{align}\n"
    "\n"
    "\n"
)
assert old_dup_expand in tex, "Duplicate expand paragraphs not found"
tex = tex.replace(old_dup_expand, "")

# ═══════════════════════════════════════════════════════════════════════════
# 5. Merge "Definitions" + "Force Balance" + "Connecting Lotter to shear"
#    into one clean appendix subsection.
#    Use index-based extraction to avoid whitespace mismatch.
# ═══════════════════════════════════════════════════════════════════════════
start_marker = r'\paragraph{Definitions}'
end_marker = '% \u2500\u2500 Supplementary figures'
i_start = tex.find(start_marker)
i_end = tex.find(end_marker)
assert i_start != -1, "Definitions paragraph not found"
assert i_end != -1, "Supplementary figures marker not found"
old_force = tex[i_start:i_end]
new_force = (
    "\\subsubsection*{Force balance and the Lotter connection}\n"
    "\n"
    "For a binary hillslope with vegetated patches (length $L_V$, "
    "roughness $n_V$, depth $h_V$) and bare patches ($L_B$, $n_B$, "
    "$h_B$), the momentum balance along the slope reads\n"
    "\\[\n"
    "\\rho (L_V h_V + L_B h_B)\\,g\\tan\\theta\n"
    "= L_V \\tau_V\\,h_V + L_B \\tau_B\\,h_B.\n"
    "\\]\n"
    "Expressing the shear stress through the Manning drag coefficient "
    "$C_d = g n^2 h^{-1/3}$ gives\n"
    "\\[\n"
    "(L_V h_V + L_B h_B)\\,g\\tan\\theta\n"
    "= L_V \\frac{C_{d,V}}{n_V^2}\\,h_V^{7/3}\n"
    "+ L_B \\frac{C_{d,B}}{n_B^2}\\,h_B^{7/3}.\n"
    "\\]\n"
    "When depth and velocity differences between patches are neglected "
    "($h_V = h_B = h$, $U_V = U_B = U$), this simplifies to an "
    "area-weighted average of the drag coefficients,\n"
    "\\[\n"
    "C_{d,\\mathrm{eff}} = \\frac{L_V C_{d,V} + L_B C_{d,B}}"
    "{L_V + L_B},\n"
    "\\]\n"
    "which, via $n = (C_d h^{1/3}/g)^{1/2}$, yields $n_{\\mathrm{eff}}"
    " \\sim h^{2/3}$ rather than the $h^{5/3}$ weighting of the Lotter "
    "equation. Allowing $h_V \\neq h_B$ restores the Lotter form, "
    "confirming that the depth covariance is the mechanism by which "
    "the equal-discharge result produces $n_e < \\langle n \\rangle$.\n"
)
tex = tex[:i_start] + new_force + tex[i_end:]

# ═══════════════════════════════════════════════════════════════════════════
# Write
# ═══════════════════════════════════════════════════════════════════════════
tex_path.write_text(tex)

old_lines = original.count('\n')
new_lines = tex.count('\n')
print(f"Patched {tex_path}")
print(f"  Lines: {old_lines} → {new_lines}  (removed {old_lines - new_lines})")
print(f"  \\includegraphics count: {tex.count('includegraphics')}")
