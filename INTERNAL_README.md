# Internal notes — Roughness Scale project

> **Not tracked in git.** For public-facing info see `README.md`.

---

## Overleaf sync

`overleaf/` is a git clone of project `663fbe37ea3d10144699322a`.

```bash
# Pull latest from Overleaf
git -C overleaf pull

# Push local edits to Overleaf
git -C overleaf add .
git -C overleaf commit -m "update figures"
git -C overleaf push
```

Files in `overleaf/`:
- `main.tex` — main manuscript body
- `equiv_rough.tex`, `spatial_averaging.tex`, `lit review.tex` — section files
- `references.bib`, `local.bib`
- `figures/` — figures included in the manuscript

---

## Project TODOs

### Manuscript

- [ ] Write **Results** section — currently a placeholder stub
- [ ] Write **Discussion** section — currently almost empty; key points:
  - Why is $n_e < \langle n \rangle$? (run-around / preferential routing)
  - Link to composite-channel literature (Flintham & Carling 1992, Djajadi 2009)
  - Compare with porous-media / canopy analogies from lit review
- [ ] Add missing citations: emulation/hybrid ML approaches (`TODO` in `main.tex` line 82)
- [ ] Decide on double-averaging vs pure spatial-averaging presentation of NS equations (`TODO` in `main.tex` line 110)
- [ ] Address **bounds questions** in Results:
  - Forward: given $f_V$, $\sigma$, $\ell$, anisotropy → predict $r_\mathrm{eq}$
  - Inverse: given calibrated $n_e$ → constrain patch statistics
  - Parallel vs series stripe configurations as theoretical bounds
- [ ] Complete literature review on overland flow roughness scaling (Kim, Thompson 2011 and others noted in `main.tex`)
- [ ] Finalize notation consistency: $r_\mathrm{eq}$ vs $n_e/\langle n \rangle$ — pick one and apply throughout

### Analysis & figures

- [ ] **Figs 1–2** not yet in registry — create and register
  - Fig 1: domain schematic / example spatial patterns
  - Fig 2: example hydrographs / IF curves showing $n_e$ fitting
- [ ] **Storm scatter plot** (roughness_scale-pattern.ipynb) — register as a numbered figure
- [ ] **Compare $n_e^\mathrm{IF}$ vs $n_e^\mathrm{Q}$** — dedicated figure showing agreement/divergence across parameter space
- [ ] Bounds figure: parallel-stripe vs series-stripe vs random at fixed $f_V$
- [ ] Update figure registry source notebook names — currently all show `roughness_scale-analysis.ipynb`; re-save from `roughness_scale-pattern.ipynb` / `roughness_scale-decomp.ipynb`

### Code / infrastructure

- [ ] `roughness_scale-pattern.ipynb` — verify all cells run clean top-to-bottom after split
- [ ] `roughness_scale-decomp.ipynb` — verify all cells run clean top-to-bottom after split
- [ ] Check whether `stats.py` functions (`fit_ols`, `eval_yhat`) are imported correctly in decomp notebook after split
- [ ] Consider extracting `plot_CF_errors` and `plot_er_ratio_4cols` to `src/plot_helpers.py`
