---
type: preregistration
project: priscus_excerpta
date: 2026-11-11
status: FROZEN before any hardware run. NOT LAUNCHED (cloud session: no Irene access). Launch needs CJ's go, from a Mac.
skill: ruby-run (read) · procedure 15 · QPU launch gate
files: ruby_mis_prep.py · layout_ruby.py · layout_ruby.json · ruby_mis_prep.json (register, preflight, emulation)
---

# Pre-registration: Priscus stemma as MIS on Ruby (no addressing)

## Question

Does Ruby, with a uniform drive and no per-atom addressing, return the maximum independent sets of the
verified Priscus family graph? These are the four most-resolved family trees the evidence allows. And in
what proportions, against emulation at Ruby geometry?

The evidence-weighted version (MWIS, de Boor's tree) is **emulation only**: Ruby has no DMM.

## Instance (frozen)

- Families (rooted shared errors, de Boor I.1 pp. 121–155, read by eye by two readers): B+M+P 45,
  M+P 8, E+B 4, E+P 3, E+M+P 3 ([[RUN_PRISCUS_VERIFIED_20261110]]).
- Conflict graph: 6 edges.
- Maximum independent sets (size 2): {B+M+P, M+P} (de Boor), {M+P, E+B}, {M+P, E+M+P}, {E+P, E+M+P}.
- Exact MWIS: {B+M+P, M+P}.

## Register and drive (frozen; `ruby_mis_prep.json`)

| | value |
|---|---|
| layout | `layout_ruby.json`: conflicts ≤ 0.80, compatible ≥ 1.30, min 0.72 (layout units); 7.2 µm/unit |
| copies | **9 copies on a ring, 45 atoms**; min spacing 5.18 µm; radius ≤ 50 µm; min inter-copy distance 13.7 µm (V ≤ 0.13 rad/µs) |
| drive | Ω: 0 → 2π·1.5 MHz → 0 (9.42 rad/µs ≤ 15.708 QPU cap); δ: −2π·4 → +2π·0.9 MHz (interpolated, 4 knots); **T = 6000 ns** |
| blockade | R_b(Ω) = 6.72 µm; conflict pairs ≤ 0.86 R_b; compatible ≥ 1.39 R_b (house margin ≥ 1.25 ✓) |
| energies | V_conflict ≥ 23.7 rad/µs > 2 δ_f (11.3); V_compatible ≤ 1.28 < δ_f/4 (1.41) |
| shots | **500** (= 4,500 copy-samples) |

All preflight asserts pass in `ruby_mis_prep.py`: Ruby QPU amplitude cap, detuning cap, duration,
5 µm spacing, 50 µm radius, blockade margin, and the final detuning against the residual tail.
The device was a RubyLike virtual copy of Pulser `AnalogDevice` (same C6 = 865,723 rad/µs·µm⁶, n = 60),
with radius 50 µm and 100 atoms. **On the Mac, rebuild against the served spec**
(`Device.from_abstract_repr(ruby_device.json)`) with `Register({...}).with_automatic_layout(device)`, and
re-run the asserts in the container before submitting (skill rules 1 and 11).

## Emulation predictions (frozen; QutipEmulator, one copy)

| | ideal | noisy (vendor nominal: p_fp 0.02, p_fn 0.10, T2* 4.5 µs) |
|---|---|---|
| MIS run: P(valid IS) | 1.000 | 0.926 |
| MIS run: P(any maximum IS) | 0.922 | 0.512 |
| MIS run: per set {B+M+P,M+P} / {M+P,E+B} / {M+P,E+M+P} / {E+P,E+M+P} | 0.225 / 0.137 / 0.447 / 0.114 | 0.124 / 0.108 / 0.166 / 0.114 |
| MWIS emulation (DMM, w_eff = 0.3 + 0.7 w/w_max, optimum unchanged): P(de Boor's tree) | 0.568 | 0.190 |

Random baseline: a uniform random 5-bit string hits one of the four maximum sets with probability
4/32 = 0.125.

## Hypotheses and pass rules (hardware, per copy, pooled over 9 copies × 500 shots)

- **H1, validity:** P(valid independent set) ≥ 0.80. PASS if the 95% Wilson lower bound is ≥ 0.80.
- **H2, MIS:** P(any maximum IS) is above the random baseline. PASS if the 95% Wilson lower bound is
  > 0.125 **and** the point estimate is ≥ 0.25.
- **H3, the space of histories:** all four maximum sets appear, each at ≥ 2% of copy-samples.

**Reported, no pass/fail:**
- per-set shares against the noisy emulation (total-variation distance);
- per-copy (ring position) variation;
- a readout-corrected secondary analysis (confusion matrix from the vendor nominal values), labelled
  as secondary.

**Not claimed:** speed, advantage, or that the hardware chooses de Boor's tree. That needs weights, and
is emulation only.

## Launch (CJ, from a Mac with `IRENE-OK`)

1. Step 0 of the `ruby-run` skill; `ruby_status.py` must not be `DOWN`; no other chain running.
2. Append ONE row to `JOB_TABLE`: register = `ruby_mis_prep.json["register_um"]` (45 atoms); the drive
   above; 6000 ns; 500 shots; label `priscus_mis_k001`. Print and assert the row inside `ccc-quantum`.
3. `./launch.sh <ROW> 500`. Score the returned file per copy with the frozen rules above. Record it in the
   QPU campaign ledger.
