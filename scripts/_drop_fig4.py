#!/usr/bin/env python3
"""Remove fig4 from main.tex, renumber fig5→fig4, fig6→fig5."""
import pathlib

tex_path = pathlib.Path(__file__).resolve().parent.parent / "overleaf" / "main.tex"
tex = tex_path.read_text()

# 1. Remove the entire fig4 figure environment
fig4_block = (
    "\n\\begin{figure}[htbp]\n"
    "  \\centering\n"
    "  \\includegraphics[width=\\textwidth]{fig4_Sf_decomp_3x3.png}\n"
    "  \\caption{As in Fig.~\\ref{fig:Sf-decomp} but for the full storm ensemble "
    "(all anisotropic patterns, $p=8$~mm/hr, $t_r=60$~min). Row~1: observed, $T_2$, "
    "and Lotter effect ratios. Rows~2--3: individual fluctuation terms, coloured by~$L_V$.}\n"
    "  \\label{fig:Sf-decomp-full}\n"
    "\\end{figure}"
)
assert fig4_block in tex, "fig4 block not found"
tex = tex.replace(fig4_block, "")

# 2. Rename fig5 → fig4, fig6 → fig5 in includegraphics
tex = tex.replace(
    "\\includegraphics[width=\\textwidth]{fig5_obs_vs_pred_re_norm_9panel.png}",
    "\\includegraphics[width=\\textwidth]{fig4_obs_vs_pred_re_norm_9panel.png}"
)
tex = tex.replace(
    "\\includegraphics[width=\\textwidth]{fig6_pred_vs_obs_CF.png}",
    "\\includegraphics[width=\\textwidth]{fig5_pred_vs_obs_CF.png}"
)

tex_path.write_text(tex)
print(f"Patched {tex_path}")
print(f"  includegraphics count: {tex.count('includegraphics')}")
print(f"  fig4_ references: {tex.count('fig4_')}")
print(f"  fig5_ references: {tex.count('fig5_')}")
print(f"  fig6_ references: {tex.count('fig6_')}")
