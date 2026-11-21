# stemma-mwis

Code and data for **C. Jurczak, *Testing stemmata against the evidence: Exact optimisation, evidence maps and Bédier's question*** (2026, preprint).

Archived release: https://doi.org/10.5281/zenodo.24405388. Technical report with the proofs and the audit trail of the checks: https://doi.org/10.5281/zenodo.24405385. Lean formalisation: https://github.com/metamorphoses-collection/stemma-lean (https://doi.org/10.5281/zenodo.24405390).

## Dataset: MGH *Auctores antiquissimi* critical apparatus in tables

**[`mgh_apparatus/`](mgh_apparatus/)** — the critical apparatus of 42 texts of the *Monumenta Germaniae Historica, Auctores antiquissimi* (1877–1919), read by machine from the page images: 122,054 entries and 150,947 readings, one CSV per text (reading, lemma, witnesses, page image, transcribed apparatus line), with tables of texts, sigla and extraction checks, and the 4,091 page transcriptions. CC BY 4.0. Also downloadable on its own from the Zenodo record as `MGH_Auctores_antiquissimi_apparatus_tables.zip`.

A stemma is reconstructed as a maximum-weight set of jointly compatible families (splits) proposed by shared readings. Fragmentary witnesses are handled exactly: missing passages become joint completion of partial splits, solved by cutting planes on top of the pairwise relaxation. The code also draws an **evidence map** of any stemma (support, conflict and bootstrap stability of each family), benchmarks language models on the same input, and prepares native encodings for neutral-atom quantum processors.

## Layout

| Folder | Contents |
|---|---|
| `code/benchmark/` | exact solver (`exact_partial.py`: overlap relaxation + CP-SAT joint-completion check + no-good cuts; deterministic by default), shared utilities (`stemma_mwis.py`), CASC skeleton scorer (`score.py`), bootstrap + parsimony protocol (`mwis_protocol.py`, `pars.R`), conflict report (`contamination_report.py`), earlier pairwise heuristic (`partial*.py`), language-model benchmark (`llm_bench.py`) and evidence checkers (`llm_check_*.py`) |
| `code/donatus/` | building a witness × unit matrix from a transcribed negative apparatus (`build_matrix.py`, spec in `TRANSCRIPTION_SPEC.md`), double-transcription checks (`proof_diff.py`, `reconcile.py`), baselines (neighbour joining, scoring against an editor's families: `baselines.py`), the repaired Camps–Cafiero routine (`pcc_fixed.R`), the Lawrence matrix builder (`lawrence_matrix.py`) |
| `code/native/` | conflict graph of a tradition, hardware-valid unit-disk layout, bootstrap over candidate histories |
| `ruby/` | neutral-atom emulation (Pulser 1.10): preflight checks against device limits, MIS and weighted MIS (local detunings), ideal and noisy, full sample distributions |
| `data/lawrence/` | BHL 4759 prologue: 7 × 75 numeric matrix, protocol outputs, native graph, layout, emulation results |
| `data/donatus_results/`, `data/heinrichi_results/`, `data/llm_results/` | derived results only (stabilities, trees with sigla, scores, evidence checks) |
| `code/historic/` | re-reading printed editions: deterministic parser for Mommsen's *Getica* apparatus (`mommsen_jordanes_extract.py`), generic extractor of shared departures from MGH-style negative apparatus (`mgh_rooted_extract.py`), rooted exact MWIS with bootstrap (`rooted_mwis.py`), scoring an extraction against page-image gold groups (`proof_score.py`), re-OCR helpers (dMGH + Tesseract: `dmgh_tesseract_page.sh`, `reocr_assemble.py`; Transkribus API: `transkribus_*.py`, credentials read from your own `~/.config/transkribus/env`), LLM-assisted extractors with verification for Lachmann's *Iwein* and Lucretius (`iwein_extract.py`, `lachmann_extract.py`) |
| `code/benchmark/network_bayes_baselines.py` | Neighbor-Net (R phangorn) and Bayesian (MrBayes 3.2.8) baselines, scored like the others |
| `data/mgh/` | five *Monumenta Germaniae Historica* editions (Mommsen's *Getica* and *Variae*, Droysen's Eutropius, Luetjohann/Leo's Sidonius, Leo's Fortunatus): extracted units and rooted results; `proof/gold_groups.json` = 200 printed groups read on 12 page images, the test used to measure extraction (recall 84%) |
| `data/mommsen_getica/` | the *Getica* with the dedicated parser: units, rooted results by Mommsen's lemma and by the majority reading, unrooted protocol |
| `data/baselines/` | Neighbor-Net and MrBayes scores and consensus trees (Heinrichi run not converged after 4×10^6 generations; consensus stable) |
| `data/llm_results/rep2`, `rep3` | second and third runs of every model |

## Sources for the printed editions

The MGH volumes (1879–1894) are in the public domain. OCR text: Internet Archive items `cuaiordanisroman00jord`, `cuaeutropibreviar00eutr`, `gaisoliiapollin00ltgoog`, `cassiodorisenato00cass`, `venantihonoricl00unkngoog` (download the `_djvu.txt` files); page images: dMGH, https://www.dmgh.de/ .

## Data not included

- **Heinrichi** (CASC 2007): download `heinrichi.num` and `correct.txt` from the challenge site, https://www.cs.helsinki.fi/u/ttonteri/casc/ , into `data/heinrichi/`. The organisers' terms apply.
- **Donatus** (BHL 2294): the matrix is derived from the apparatus of P. Licciardello's edition (SISMEL, 2018) and is not redistributed. `code/donatus/build_matrix.py` rebuilds it from a transcription made by the user, following `TRANSCRIPTION_SPEC.md`.
- **Language-model prompts** embed the Heinrichi matrix and are regenerated by `llm_bench.py` from the downloaded data.

## Reproducing the main results

```
pip install ortools numpy scipy networkx      # Python 3.12; R with phangorn for parsimony
python code/benchmark/exact_partial.py data/heinrichi/heinrichi.num --truth data/heinrichi/correct.txt --out out/h
python code/benchmark/mwis_protocol.py data/lawrence/lawrence_bhl4759.num --out out/lawrence --B 200 --thetas 0.5,0.7,0.9
python code/native/lawrence_graph.py data/lawrence/lawrence_bhl4759.num out/graph.json 3
python code/native/bootstrap_histories.py data/lawrence/lawrence_bhl4759.num out/graph.json 20000
```

The exact solver is deterministic by default (one CP-SAT worker, seed 0); set `MWIS_WORKERS` and `MWIS_SEED` to explore co-optimal solutions.

## Pre-registered test of the editors' top divisions (MGH, Sec. "A pre-registered test across the series")

- `data/bedier/docs/` — the frozen protocol with all amendments (`PREREG_BEDIER_20261113.md`), results with the
  exploratory analyses (`RESULTS_BEDIER_20261113.md`), the survey of *Scriptores rerum Merovingicarum*
  (`SURVEY_BEDIER_SRM_20261113.md`) and the survey of later editions (`LATER_EDITIONS_SURVEY_20261113.md`).
- `data/bedier/tops/`, `data/bedier_srm/tops/` — editors' stated top divisions, recorded blind from the prefaces
  (MGH *Auctores antiquissimi* 1–15, 71 texts; *Scriptores rerum Merovingicarum* 1–7, 168 texts).
- `data/bedier/ocr_llm/` — vision-model transcriptions (google/gemini-4-argon, temperature 0, prompt in
  `code/historic/bedier_llm_read.py`) of the 4,091 dMGH pages of the located texts, one file per page
  (`<volume>_<leaf>.llm.txt`); `USAGE.tsv` logs tokens and cost. Five blank pages are marked as such.
  These are machine transcriptions of public-domain 19th-century editions, checked on sampled pages only.
- `data/bedier/gold/` — blind gold transcriptions of the proof, development and fresh proof pages (AI agents reading
  the page images; not checked by a human), with the page draws (`leaves*.json`, `a10_draws.json`).
- `data/bedier/runs/` — extracted units per text (parser v4); `screen_a10/`, `classification_a10/` — the frozen
  screen and the single classification run; `proof_scores*.json`; `null_exploratory_primary.json`,
  `why_exploratory.json` — exploratory analyses.
- `code/historic/` — `bedier_a6.py` (locate texts, extract), `bedier_llm_read.py`, `mgh_rooted_extract.py` (parser v4),
  `bedier_proof.py`, `bedier_top.py` (frozen scorer), `bedier_null.py`, `bedier_why.py`.

The Apple Vision / Tesseract OCR of all 9,385 pages of MGH AA 1–15 (four variants per page, amendment A6) is attached
to the GitHub release as `mgh_aa_apple_tesseract_ocr_4variants.zip`. The dMGH page images are not redistributed;
they are at https://www.dmgh.de/ (`<volume>/img/300/<volume>_<leaf>.jpg`).

Reproduce the classification:
```
python code/historic/bedier_top.py data/bedier/runs/<text>/units.json --witnesses <sigla> --top '<A|B>' --out out.json --B 200
```
(sigla and top per text in `data/bedier/settings.json`).

## Bédier's own case: the *Lai de l'Ombre* (Sec. "Bédier's own case")

- `data/lai_ombre/` — vision-model transcriptions of the variant notes of Bédier's second edition (1913, Internet Archive
  `lelaidelombre00jeanuoft`; `vision/`, one file per page), the 814 extracted notes (`units.json`), the blind proof page
  (`proof_draw.json`, `proof_gold_n103.json`, `proof_score_n103.json`: 42 of 43 readings), and the three historical
  classifications scored against the readings (`class_bedier1890.json`, `class_paris1890.json`); `NOTE_LAI_OMBRE.md`
  summarises the analysis. Code: `code/historic/lai_read.py`, `lai_extract.py` (OpenRouter key read from your own environment).

## Contamination: which witness carries the conflict (exploratory)

`code/historic/contam_witness.py` deletes one witness at a time, re-solves the rooted optimum and compares the drop in
conflict weight with the witness's share of the evidence, with a unit bootstrap. Results for Ennodius and the *Getica*:
`data/bedier/contam_exploratory.json`.

## Trees

`data/trees/` — the stable families of every classified text as rooted Newick trees (labels: family, support, stability),
and the unrooted tree of the *Lai de l'Ombre*. Built by `code/historic/families_to_newick.py`.

## Neutral-atom runs (Sec. "Stemmata on neutral-atom quantum processors")

- `data/ruby/priscus/` — the frozen protocol and its amendment, the register and drive (`ruby_mis_prep.json`), the
  scorer (`score_ruby_mis.py`) and the four hardware result files of the run on Pasqal's Ruby (rows 728–731,
  3,500 shots); score them with `python score_ruby_mis.py ruby_mis_prep.json ruby_diag_job_72{8,9}.json ruby_diag_job_73{0,1}.json`.
- `data/ruby/lawrence/`, `data/ruby/lai/` — frozen protocols, registers, emulated predictions, scorers and the hardware
  result files for the Lawrence core (row 732, 5 families, 7 copies, 35 atoms, 1,000 shots) and the *Lai de l'Ombre*
  (row 733, 7 families, 8 copies, 56 atoms, 1,000 shots), with their scores. Score them with
  `python score_ruby_mis_lawrence.py lawrence_mis_prep.json ruby_diag_job_732.json` and
  `python score_ruby_mis_lai.py lai_mis_prep.json ruby_diag_job_733.json`. `shot_bootstrap_lawrence.py` is a post-hoc
  sensitivity check of the Lawrence sampler share that resamples whole shots (seven copies together); it does not change
  the frozen analysis.
- `ruby/` — register construction and emulation (`ruby_mis_prep.py`, `emu_dist.py`, `layout_ruby.py`,
  `lai_layout_compact.py`, `lai_tile.py`).

## Licence

Code: MIT (see `LICENSE`). Data files in `data/`: CC BY 4.0.

## Mathematical checks, Fable runs and supplement

- `identifiability/`: scripts and outputs of the machine checks C1–C25 and the simulations behind the
  technical report (https://doi.org/10.5281/zenodo.24405385) and the mathematical paper (see
  `identifiability/README.md`).
- `data/llm_results/fable_runs/`: the three Heinrichi runs and three Donatus runs of Claude Fable 5.2
  reported in the paper (every answer kept in full).
- `papers/supplement_testing_stemmata.pdf`: Supplementary material S1–S10 of *Testing stemmata
  against the evidence*.
- Lean formalisation: https://github.com/metamorphoses-collection/stemma-lean.
