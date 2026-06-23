# MWQTM — JASA-structured manuscript

Event-aligned conditional **quantile** trajectories of a biomarker before disease
onset under left truncation, right censoring, competing death, and irregular
visits. This repository is a **JASA-style restructure** of the original MWQTM
draft (kept intact at the companion `MWQTM` repository): a lean main paper
(*model → identification → orthogonal inference → concise simulations →
application*) with all proofs, algorithms, audit/sensitivity studies, and the
Sun et al. comparison moved to a supplement.

## Layout

```
paper/
  paper.tex            main manuscript (~9 pp): Intro, Model, Estimation,
                       Theory, Simulations, Application, Discussion
  supplement.tex       supplement (~43 pp): A proofs · B algorithms ·
                       C additional simulations · D application details
  proof_body.tex, bahadur_body.tex, censoring_orthogonality_body.tex
                       imported proof bodies (preamble-stripped)
  table_*.tex          audit/sensitivity/comparison tables (all -> supplement)
  figs/                cmp_*.pdf (Sun comparison), app_panels.pdf (application)
  mwqtm_references.bib
results/               figures reused from the simulation pipeline
R/                     estimator + simulation + comparison code (reproducibility)
```

## Build

```
cd paper
pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
pdflatex supplement && pdflatex supplement      # xr cross-refs the main paper
```

## Notes

- **One formal inferential route in the main text:** the cross-fitted one-step
  estimator (B-OS-CF, Theorem 1) carries the √n asymptotic-normality result; the
  joint-sieve A-CV estimator is the practical point estimator (subject bootstrap),
  shown in simulation to agree with B-OS-CF.
- **Sun et al. comparison** (Section 5.4 / Supplement C): a point-estimation
  comparison against the *actual* estimator of Sun et al. (2025); no inference
  claim is made for the recovered mean.
- **Application** (Section 6): an *illustrative analysis on a synthetic
  ADNI-matched cohort* — clearly labeled, not real data. The substantive ADNI
  analysis is in preparation and will replace it.
