#!/usr/bin/env python3
"""Insert channel-equations table into overleaf/main.tex, inside the Bridge subsection."""
import pathlib

tex_path = pathlib.Path(__file__).resolve().parent.parent / "overleaf" / "main.tex"
tex = tex_path.read_text()

# ── 1. Insert the table after the Bridge introductory paragraphs ────────
# Place the table right before the first \begin{figure} in the Bridge section.
anchor = (
    "Furthermore, in our case the averaging is along the line of flow, "
    "whereas the composite-channel problem concerns the cross-sectional "
    "distribution of roughness (e.g., bed versus side slopes with "
    "different materials).\n"
    "\n"
    "\\begin{figure}[htbp]"
)

table_tex = r"""Furthermore, in our case the averaging is along the line of flow, whereas the composite-channel problem concerns the cross-sectional distribution of roughness (e.g., bed versus side slopes with different materials).

Table~\ref{tab:channel-eqs} summarises the nine composite-channel formulas tested. In each formula, $n_i$ and $h_i$ denote the local Manning coefficient and flow depth in sub-section~$i$, angle brackets $\langle\cdot\rangle$ denote the area-weighted spatial average (which, for uniform-width overland flow, equals the arithmetic mean over cells), and $\bar{h} = \langle h_i \rangle$.

\begin{table}[htbp]
\centering
\small
\caption{Composite-channel equations for equivalent Manning's $n$. All formulas reduce to area-weighted averages for uniform-width overland flow. References marked $^\dagger$ are cited via \citet{djajadi2009comparative}.}
\label{tab:channel-eqs}
\begin{tabular}{c l l l}
\hline
\# & Method & Formula & Reference \\
\hline
1 & Lotter          & $\displaystyle n_c = \frac{\bar{h}^{5/3}}{\left\langle \dfrac{h_i^{5/3}}{n_i} \right\rangle}$                   & \citet{lotter1933considerations} \\[14pt]
2 & Pavlovskii      & $\displaystyle n_c = \sqrt{\langle n_i^{2} \rangle}$                                                            & Pavlovskii (1931)$^\dagger$ \\[10pt]
3 & Horton--Einstein & $\displaystyle n_c = \langle n_i^{3/2} \rangle^{2/3}$                                                          & Horton (1933)$^\dagger$; Einstein (1934)$^\dagger$ \\[10pt]
4 & Cox (arithmetic) & $\displaystyle n_c = \langle n_i \rangle$                                                                      & \citet{cox1973effective} \\[10pt]
5 & Krishnamurthy--Christensen & $\displaystyle n_c = \exp\!\left(\frac{\langle h_i^{3/2} \ln n_i \rangle}{\langle h_i^{3/2} \rangle}\right)$ & K.--C.\ (1972)$^\dagger$ \\[14pt]
6 & Colebatch        & $\displaystyle n_c = \left(\frac{\langle h_i\, n_i^{3/2} \rangle}{\bar{h}}\right)^{\!2/3}$                     & Colebatch (1941)$^\dagger$ \\[14pt]
7 & Yen $R^{1/6}$    & $\displaystyle n_c = \bar{h}^{1/6}\,\left\langle \frac{n_i}{h_i^{1/6}} \right\rangle$                          & \citet{yen1992channel} \\[14pt]
8 & Yen $R^{1/3}$    & $\displaystyle n_c = \frac{\langle h_i^{1/3}\, n_i \rangle}{\bar{h}^{1/3}}$                                    & \citet{yen1992channel} \\[10pt]
9 & Felkel           & $\displaystyle n_c = \frac{1}{\langle n_i^{-1} \rangle}$                                                       & Felkel (1960)$^\dagger$ \\[10pt]
\hline
\end{tabular}
\end{table}

\begin{figure}[htbp]"""

assert anchor in tex, f"Anchor not found in main.tex"
tex = tex.replace(anchor, table_tex, 1)

tex_path.write_text(tex)
print(f"Patched {tex_path}")
print(f"  \\begin{{table}} count: {tex.count(chr(92) + 'begin{table}')}")
print(f"  Total lines: {len(tex.splitlines())}")
