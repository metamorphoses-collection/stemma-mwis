---
type: preregistration
project: priscus_excerpta
date: 2026-11-14
status: FROZEN 2026-11-14, before any hardware run (CJ, 14 Nov: "if we can do as good as priscus then we should have le lai
  on ruby, prepare and run after the current runs (active)"). Irene JOB_TABLE row 733 `lai_mis_k001`, 1000 shots, queued
  afterany the last active chain (one chain at a time). Scorer score_ruby_mis_lai.py self-tested on synthetic results drawn
  from both emulations (noiseless and noisy: H1-H3 PASS).
files: data/lai_layout_ruby.json (graph) · data/lai_layout_ruby_s11.json (single-copy layout, lai_layout_compact.py) ·
  data/lai_ruby_prep/ruby_mis_prep.json (ruby_mis_prep.py: preflight + emulation) · ruby/lai_mis_prep.json (8-copy register,
  lai_tile.py) · data/lai_ruby_dist/emu_dist.json (full emulated distributions at the frozen sweep, emu_dist.py)
---

# Pre-registration: Bédier's *Lai de l'Ombre* on Ruby

## Question
Does Ruby, with the uniform drive of the Priscus rows 728–731, return the maximum independent sets of the *Lai de l'Ombre*
family graph as well as it did for Priscus? This is the text of Bédier's 1928 objection (paper Sec. "Bédier's own case").
Uniform drive gives the MIS, not the MWIS: the four maximum sets are tied, and the device is not expected to prefer the
edition's tree among them.

## Instance (frozen with this file)
- Units: data/historic/lai_ombre/units.json (Bédier 1913 apparatus, proof page passed). Split weights on the seven
  manuscripts A–G, shared departures of 2–5 manuscripts counted as splits.
- **Support ≥ 11 units**: seven families CG 176, DE 43, AB 38, DEF 24, EF 23, CEG 13, DF 13. (At ≥ 10, ABF 10 enters and
  makes an induced C4 that has no unit-disk embedding with these margins; ≥ 11 is the smallest threshold that embeds.)
- Conflicts (incompatible splits), 6 edges: DE–EF, DE–CEG, DE–DF, DEF–CEG, EF–CEG, EF–DF. AB and CG conflict with nothing.
- Maximum independent sets (size 4, tied): {CG,DE,AB,DEF}, {CG,AB,DEF,EF}, {CG,AB,DEF,DF}, {CG,AB,CEG,DF}.
- Unique MWIS {CG,DE,AB,DEF}, weight 281 (next 261): the unrooted tree ((A,B),(C,G),((D,E),F)) on which Bédier 1890 and
  Quentin agree.

## Register and drive (frozen; drive identical to rows 728–731)
- One copy: 7 atoms, conflicting pairs ≤ 0.847 R_b (5.69 µm), compatible pairs ≥ 1.404 R_b (9.43 µm), min spacing
  5.26 µm, R_b = 6.72 µm at Ω = 9.42 rad/µs (C6 865,723 rad/µs·µm⁶). V_conflict ≥ 25.6, V_compatible ≤ 1.23 rad/µs.
- **8 copies × 7 atoms = 56 atoms**, hexagonal tiling, inter-copy ≥ 16.4 µm (tail ≤ 0.05 rad/µs), radius ≤ 49.4 µm.
- Global drive: Ω 0 → 9.4248 → 9.4248 → 0 rad/µs; δ −25.1327 → 5.6549 rad/µs; T = 6000 ns.
- **One row, 1000 shots** (8000 copy-samples). Register from the coordinates (`with_automatic_layout`), asserted inside
  the container against the served spec, and the drive asserted equal to row 728, before submission.

## Emulated predictions (computed before any hardware data; single copy, 4000 samples, emu_dist.json)
| | noiseless | noisy (vendor nominal: readout 0.02/0.10, T2* 4.5 µs) |
|---|---|---|
| valid independent sets | 1.000 | 0.933 |
| any maximum set | 0.738 | 0.213 |
| {CG,DE,AB,DEF} (edition tree) | 0.280 | 0.067 |
| {CG,AB,DEF,EF} | 0.046 | 0.035 |
| {CG,AB,DEF,DF} | 0.309 | 0.064 |
| {CG,AB,CEG,DF} | 0.103 | 0.048 |
| random configuration | any maximum set 4/128 = 0.031 | |
The noisy prediction is lower than for Priscus (0.51) because each maximum set has four excited atoms, so readout loss alone
costs about a third. On Priscus the hardware did better than the noisy emulation (0.62 vs 0.51, ratio 1.22).

## Primary criteria (per copy-sample, pooled over 8 copies, Wilson 95% intervals); same form as Priscus
- **H1** valid independent sets: Wilson lower bound ≥ 0.80.
- **H2** maximum sets: Wilson lower bound of P(any MIS) > 4/128, and point P(any MIS) ≥ 0.10 (half the noisy prediction).
- **H3** the space of histories: each of the four maximum sets ≥ 0.02 of copy-samples.
"As good as Priscus" is read as all three passing. Not a criterion: the hardware/noisy ratio, compared with Priscus's 1.22.

## Reported in addition (not criteria)
P(edition tree); hardware/noisy ratio of P(any MIS); total-variation distance of the four MIS shares to each emulation;
excitation of the isolated atoms AB and CG (noisy emulation 0.66 each); per-copy spread; top outcomes. The MWIS step is not
run (Ruby has no local addressing). With the Priscus weight compression (w_eff = 0.3 + 0.7 w/w_max) the weighted emulation
does not raise the edition tree here (noiseless 0.28 → 0.13), because CG 176 flattens the other weights; it is reported as
such and not used.

## Analysis plan
Score with `ruby/score_ruby_mis_lai.py ruby/lai_mis_prep.json <ruby_diag_job_733.json>`. The result file is the only truth:
audit job_idx, block and batch_id against the filename. Record in RUN_RUBY_MIS_LAI_<date>.md, the QPU ledger and NOW.md.

## Do not
Change the drive, threshold or copy count after seeing data; add rows without an amendment written before their data;
resubmit a live batch.
