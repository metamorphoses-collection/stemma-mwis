---
type: preregistration
project: priscus_excerpta
date: 2026-11-13
status: FROZEN 2026-11-13, before any hardware run (CJ, 13 Nov: "launch lawrence row on ruby"). Irene JOB_TABLE row 732
  `lawrence_mis_k001`, 1000 shots, queued afterany the QGF-H chain (one chain at a time). Scorer score_ruby_mis_lawrence.py
  self-tested on synthetic results from the emulations (noisy: H1-H3 PASS, S mismatch, share 0.674; noiseless: share 0.956).
files: data/lawrence_layout_ruby_s3.json (layout) · data/lawrence_emulation_ruby_s3.json (register, preflight) ·
  data/lawrence_emu_dist_T6000.json (emulated distributions at the frozen sweep, computed 13 Nov with ruby/emu_dist.py)
---

# Pre-registration (draft): Lawrence prologue core on Ruby, and the "sampler of histories" test

## Questions
1. **Second tradition on hardware.** Does Ruby, with a uniform drive and no addressing, return the maximum independent
   set of the native Lawrence core, the edition's stemma? Here the largest and the best-supported set coincide, so no
   per-atom detunings are needed. It would be the second real tradition on the processor, after Priscus (rows 728–731).
2. **Sampler of histories.** Two histories compete on this core. Does the device split its samples between them in the
   proportion the bootstrap gives (0.795)? The paper says so far that one emulated match is "not yet evidence" (Sec. VI C).
   This run is the pre-registered test.

## Instance (frozen with this file)
- Lawrence prologue *In divinis laudibus*, native core at support threshold s = 3 (paper Sec. VI C):
  five families, with weights in units: B+F+Mi+P+R 11, B+R 4, P+R 4, B+F+P+R 3, B+F+R 3.
- Conflicts: P+R × B+R, P+R × B+F+R (2 edges).
- **Edition's history** E = {B+F+Mi+P+R, B+R, B+F+P+R, B+F+R}: the unique maximum independent set, and also the
  maximum-weight set (21).
- **Rival history** R = {B+F+Mi+P+R, P+R, B+F+P+R} (weight 18): P+R in place of B+R and B+F+R.
- Bootstrap (75 units, 2×10⁴ replicates, conditional on the five native families): E optimal in **0.795** of replicates.

## Register and drive (frozen; identical drive to Priscus rows 728–731)
- 7 copies × 5 atoms = 35 atoms (register_um in data/lawrence_emulation_ruby_s3.json). Conflicting pairs ≤ 0.917 R_b,
  compatible pairs ≥ 1.418 R_b, inter-copy ≥ 15.9 µm. R_b = 6.72 µm at Ω = 9.42 rad/µs.
- Global drive: Ω 0 → 9.4248 → 9.4248 → 0 rad/µs; δ −25.1327 → 5.6549 rad/µs; T = 6000 ns, as in rows 728–731.
- **One row, 1000 shots** (7000 copy-samples). Same builder as Priscus, register from the coordinates
  (`with_automatic_layout`), asserted inside the container against the served spec before submission.

## Emulated predictions at exactly this setting (computed before any hardware data; 4000 samples each)
| | noiseless | noisy (vendor nominal: readout 0.02/0.10, T2* 4.5 µs) |
|---|---|---|
| valid independent sets | 1.000 | 0.952 |
| P(E) (edition's stemma = MIS) | 0.685 | 0.140 |
| P(R) (rival) | 0.032 | 0.071 |
| share s = P(E) / (P(E) + P(R)) | 0.956 | 0.664 |
| random configuration | P(E) = 1/32 = 0.031 | |
On Priscus the hardware did better than the noisy emulation (maximum sets 0.62 vs 0.51). The outcome is therefore open
between the two columns.

## Primary criteria (hardware fidelity), scored per copy-sample and pooled, Wilson 95% intervals
- **H1** valid independent sets: Wilson lower bound ≥ 0.80.
- **H2** the edition's stemma: Wilson lower bound of P(E) > 1/32, and point P(E) ≥ 0.07 (half the noisy prediction).
- **H3** ordering: P(E) > P(R), with the 95% interval of the difference above 0.

## Sampler test (stated now, with its expected outcome)
- **S** share s = P(E) / (P(E) + P(R)), with a Wilson 95% interval on the E-or-R subsample.
  **Match** if the interval lies within [0.745, 0.845] (bootstrap 0.795 ± 0.05). **Mismatch** if the interval excludes
  that band. **Undetermined** otherwise.
- **Expected under the nominal noise model:** mismatch (predicted 0.66). A match would be reported as a match on one
  instance at one setting. It is not evidence that the device samples stemmatic uncertainty in general; that claim
  would need several instances with fixed settings. A mismatch is reported as such, with s compared with both emulated
  predictions.

## Reported in addition (not criteria)
Full outcome distribution; per-copy spread; total-variation distance to each emulated distribution; drift against the
Priscus rows (same drive), if relevant.

## Analysis plan
Score with `ruby/score_ruby_mis_lawrence.py ruby/lawrence_mis_prep.json <result.json>` (written and self-tested before
submission). The result file is the only truth: audit job_idx, block and batch_id against the filename.
Record the run in RUN_RUBY_MIS_LAWRENCE_<date>.md, the QPU ledger and NOW.md.

## Do not
Change the drive or the copy count after seeing data; add rows without an amendment written before their data;
resubmit a live batch.
