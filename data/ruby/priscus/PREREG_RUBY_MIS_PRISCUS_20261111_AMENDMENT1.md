---
type: preregistration-amendment
project: priscus_excerpta
amends: PREREG_RUBY_MIS_PRISCUS_20261111.md (frozen, untouched)
date: 2026-11-11 (13:05 CET, Studio)
status: written and submitted BEFORE any Priscus hardware data existed (row 728 still queued, Slurm 6218238 pending)
---

# Amendment 1: three more passes, same sequence

CJ, 11 Nov: "a bit more on priscus so it has good stats."

Rows 729, 730, 731 (`priscus_mis_k002`, `k003`, `k004`) are the identical row 728 sequence (same builder, same
45-atom register `register_um`, same drive, 6000 ns; `assert_priscus2.py` asserted the abstract representation equal
to row 728 against the served HPCQS1 spec, ASSERT_OK 732), 1000 shots each, new indices so no result file is
overwritten. Slurm 6218330, afterany 6218238.

Analysis: the four result files are concatenated (3,500 shots = 31,500 copy-samples) and scored once under the
frozen H1–H3 rules of the pre-registration; the rules, the register and the drive are unchanged. Row 728 alone
is also reported (the pre-registered 500-shot design) so that the pooled verdict can be checked against it; a
pass-to-fail or fail-to-pass change between the two is reported, not resolved by choosing.

Reported in addition: per-pass drift (the four shares per pass, total-variation distance between passes), which
the single row could not give.
