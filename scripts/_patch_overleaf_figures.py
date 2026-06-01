#!/usr/bin/env python3
"""Patch overleaf/main.tex to add figure environments and fix SI numbering."""
import pathlib

tex_path = pathlib.Path(__file__).resolve().parent.parent / "overleaf" / "main.tex"
tex = tex_path.read_text()

# ── 1. Add \graphicspath and remove premature SI numbering ──────────────
old_preamble = (
    "\\usepackage{rotating}\n"
    "\n"
    "% ----- Bibliography: author\u2013year -----\n"
    "\\bibliographystyle{elsarticle-harv}\n"
    "\\biboptions{authoryear,sort&compress}\n"
    "\n"
    "% SI numbering\n"
    "\\renewcommand\\thefigure{S\\arabic{figure}}\n"
    "\\renewcommand\\thetable{S\\arabic{table}}\n"
    "\\setcounter{figure}{0}\n"
    "\\setcounter{table}{0}"
)
new_preamble = (
    "\\usepackage{rotating}\n"
    "\\usepackage{graphicx}\n"
    "\\graphicspath{{figures/}}\n"
    "\n"
    "% ----- Bibliography: author\u2013year -----\n"
    "\\bibliographystyle{elsarticle-harv}\n"
    "\\biboptions{authoryear,sort&compress}"
)
assert old_preamble in tex, "Preamble pattern not found"
tex = tex.replace(old_preamble, new_preamble)

# ── 2. Insert 6 main figures after the one-line Results stub ────────────
old_results = (
    "\\section{Results}\n"
    "\n"
    "Spatially averaged roughness components predict greater roughness "
    "than equivalent roughness obtained through calibration.\n"
    "\n"
    "\\section{Discussion}"
)
new_results = r"""\section{Results}

Spatially averaged roughness components predict greater roughness than equivalent roughness obtained through calibration.

\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{fig1_conceptual_sketch.png}
  \caption{Conceptual sketch of the equivalent-roughness problem. (a)~Three-dimensional view of a patchy hillslope with vegetation fraction $f_V=0.4$, $n_{\mathrm{veg}}=0.20$, $n_{\mathrm{bare}}=0.05$. (b)~Stack of uniform-roughness simulations (library), colour-coded from high~$n$ (dark green) to low~$n$ (light green), with the best-match layer highlighted. (c)~Outlet hydrograph: heterogeneous simulation (blue), uniform-roughness library members (green shades), and the matched equivalent roughness $n_e$ (green dashed).}
  \label{fig:conceptual}
\end{figure}

\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{fig2_effect_ratio_fV_2x2.png}
  \caption{Effect ratio $n_e/\langle n\rangle$ versus vegetation fraction~$f_V$. Top row: least intense storm (left) and most intense storm (right), coloured by patch lengthscale~$L_V$. Dashed line marks unity. Bottom row: outlet hydrographs at $f_V\approx 0.5$ for both storms, coloured by~$L_V$. The ratio is below unity for all cases, reaches its minimum near $f_V\approx 0.5$, and increases toward unity with increasing patch scale.}
  \label{fig:effect-ratio}
\end{figure}

\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{fig3_Sf_decomp_3x3.png}
  \caption{Spatial decomposition of the friction slope~$S_f$ for a fixed storm ($p=8$~mm/hr, $t_r=60$~min, anisotropic patterns). Row~1: observed $n_e/\langle n\rangle$, second-order prediction ($T_2$), and Lotter prediction versus~$f_V$, coloured by~$L_V$. Rows~2--3: the six second-order fluctuation terms ($\langle r^2\rangle$, $\langle\upsilon^2\rangle$, $4\langle r\upsilon\rangle$, $-\tfrac{8}{3}\langle r\eta\rangle$, $-\tfrac{8}{3}\langle\eta\upsilon\rangle$, $\tfrac{14}{9}\langle\eta^2\rangle$) versus~$f_V$. Row~1 shares a common $y$-axis; rows~2--3 share a separate common $y$-axis to allow direct comparison of term magnitudes.}
  \label{fig:Sf-decomp}
\end{figure}

\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{fig4_Sf_decomp_3x3.png}
  \caption{As in Fig.~\ref{fig:Sf-decomp} but for the full storm ensemble (all anisotropic patterns, $p=8$~mm/hr, $t_r=60$~min). Row~1: observed, $T_2$, and Lotter effect ratios. Rows~2--3: individual fluctuation terms, coloured by~$L_V$.}
  \label{fig:Sf-decomp-full}
\end{figure}

\section{Discussion}"""
assert old_results in tex, "Results section pattern not found"
tex = tex.replace(old_results, new_results)

# ── 3. Insert fig5 & fig6 after the Bridge-to-composite subsection text ─
old_bridge_end = (
    "Furthermore, in our case the averaging is along the line of flow, "
    "whereas the composite-channel problem concerns the cross-sectional "
    "distribution of roughness (e.g., bed versus side slopes with "
    "different materials).\n"
    "\n"
    "\\subsection{The force-balance perspective"
)
new_bridge_end = r"""Furthermore, in our case the averaging is along the line of flow, whereas the composite-channel problem concerns the cross-sectional distribution of roughness (e.g., bed versus side slopes with different materials).

\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{fig5_obs_vs_pred_re_norm_9panel.png}
  \caption{Normalised comparison of simulated $n_e/\langle n\rangle$ versus predicted $n_c/\langle n\rangle$ for nine composite-channel equations (Lotter, Pavlovskii, Horton--Einstein, Cox, Krishnamurthy--Christensen, Colebatch, Yen~$R^{1/6}$, Yen~$R^{1/3}$, and Felkel), coloured by patch lengthscale~$\sigma$. Lotter performs best ($R^2 = 0.79$, RMSE${} = 0.053$); Yen~$R^{1/6}$ is weakest ($R^2 = 0.04$).}
  \label{fig:channel-eqs}
\end{figure}

\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{fig6_pred_vs_obs_CF.png}
  \caption{Predicted versus observed $n_e/\langle n\rangle$ for six prediction methods. Top row: simulation-derived covariances (slope-integrated leading term, $T_0$ actual, full second-order actual). Bottom row: correction-factor approximations. Markers coloured by~$\sigma$. Full second-order with actual covariances achieves the best fit ($R^2 = 0.77$); the CF-based full second-order is weaker ($R^2 = 0.50$), indicating that analytical CF estimates introduce error at extreme patch scales.}
  \label{fig:CF-pred}
\end{figure}

\subsection{The force-balance perspective"""
assert old_bridge_end in tex, "Bridge subsection end pattern not found"
tex = tex.replace(old_bridge_end, new_bridge_end)

# ── 4. Insert SI figures + numbering before \bibliography ───────────────
old_bib = "\\bibliography{local.bib, references.bib}\n\n\n\n\\end{document}"
new_bib = r"""% ── Supplementary figures (SI numbering) ──────────────────────────────────
\renewcommand\thefigure{S\arabic{figure}}
\setcounter{figure}{0}

\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{SI1_obs_vs_pred_re_abs_9panel.png}
  \caption{Absolute comparison of observed $n_e$ versus predicted $n_c$ for nine composite-channel equations, coloured by~$\sigma$. Lotter performs best ($R^2 = 0.98$, RMSE${} = 0.007$); Colebatch is weakest ($R^2 = 0.81$).}
  \label{fig:SI-channel-abs}
\end{figure}

\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{SI2_channel_eqs_bias_hist.png}
  \caption{Relative bias histograms ($n_c/n_e - 1$) for each composite-channel method. Positive values indicate over-prediction; negative values indicate under-prediction. Red dashed lines mark the median bias.}
  \label{fig:SI-bias-hist}
\end{figure}

\bibliography{local.bib, references.bib}



\end{document}"""
assert old_bib in tex, "Bibliography/end pattern not found"
tex = tex.replace(old_bib, new_bib)

# ── Write ────────────────────────────────────────────────────────────────
tex_path.write_text(tex)
print(f"Patched {tex_path}")
print(f"  includegraphics count: {tex.count('includegraphics')}")
print(f"  Total lines: {len(tex.splitlines())}")
