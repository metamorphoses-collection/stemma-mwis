# Checks and simulations for the mathematical results

Scripts and outputs behind the technical report *Best-supported stemmata: compatibility results for
textual traditions, with proofs and checks* (https://doi.org/10.5281/zenodo.24405385) and the paper *Coincident scribal error and the Bédier
effect*. The report's audit trail (Appendix A, checks C1–C25) names each script and output here.

- `code/run_checks.sh` runs every check and reports each exit status; `code/verify_claims.py` runs
  C1–C14 and writes `data/verify_claims.json`; C15–C25 are separate scripts named in the audit trail.
- `code/sim_identifiability.py` (quartet and ten-witness simulations, output in `data/run2/`),
  `code/bedier_on_birthdeath.py` and `code/bedier_calibrated.py` (birth–death simulations, outputs
  `data/bedier_birthdeath*.json`, `data/bedier_calibrated*.json`).
- Logs of the last full run: `data/check_C*.log`.

`graph_stats.py` and `retention_test.py` read collation matrices of real traditions; the Donatus
matrix derives from a copyrighted apparatus and is not redistributed (see the main README for how
to rebuild it). The formal proofs in Lean are in
https://github.com/metamorphoses-collection/stemma-lean.
