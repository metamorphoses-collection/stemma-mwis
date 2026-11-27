# Checks and simulations for the mathematical results

Scripts and outputs behind the technical report *Testing stemmata against the evidence: Exact optimisation, evidence maps and Bédier's question* (technical report) (https://doi.org/10.5281/zenodo.24405385) and the paper *Coincident scribal error and the Bédier
effect*. The report's audit trail (Appendix A) names each script and output here.

- `code/run_checks.sh` runs the final audit, `code/audit.py`: for each result of the report, keyed by its
  LaTeX label and in report order (plus a block for the mathematical paper's results, keys `M:prop4.4` etc.),
  it runs the check(s) that test it and records what is tested, the instances, the outcome with its key
  numbers, and the script. Outputs `data/final_audit_YYYYMMDD.json` and `data/final_audit_table.tex` (a
  longtable with rows keyed by `\ref{label}`); exits non-zero on any failure. The existing scripts are the
  machinery and run unchanged (`code/verify_claims.py` writes `data/verify_claims.json` for the
  combinatorial checks formerly numbered C1–C14; C15–C26 are separate scripts); `code/audit_new_checks.py`
  holds the checks added for results that had none (also runnable on its own).
- `code/sim_identifiability.py` (quartet and ten-witness simulations, output in `data/run2/`),
  `code/bedier_on_birthdeath.py` and `code/bedier_calibrated.py` (birth–death simulations, outputs
  `data/bedier_birthdeath*.json`, `data/bedier_calibrated*.json`).
- Logs of the last full run: `data/check_*.log`.

`graph_stats.py` and `retention_test.py` read collation matrices of real traditions; the Donatus
matrix derives from a copyrighted apparatus and is not redistributed (see the main README for how
to rebuild it). The formal proofs in Lean are in
https://github.com/metamorphoses-collection/stemma-lean.
