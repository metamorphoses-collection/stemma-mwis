# Analysis of the three Ruby runs

Classical baselines at equal shots, shot-clustered uncertainties, per-copy shares, the Lawrence optimum/rival share,
the weighted Priscus emulation and the end-of-sweep energy with all residual interactions.

The scripts expect the layout of the run folder: the prep registers (`ruby_mis_prep.json`, `lawrence_mis_prep.json`,
`lai_mis_prep.json`) beside them and the hardware results in `results/` (`ruby_diag_job_728..733.json`; copies in
`../priscus`, `../lawrence`, `../lai`). Outputs: `baselines.json`, `uncertainty.json`, `emulation_and_tails.json`.
