# Results: does each MGH editor's own apparatus support his top division? (Bédier test)
Protocol: PREREG_BEDIER_20261113.md (frozen c9089723), amendments A1–A10, all written before the data they govern.
Classification run once, 2026-11-13, on the final set (data/historic/bedier/classification_a10/; scorer bedier_top.py,
exact rooted MWIS, B = 200, stable ≥ 0.90, contradicting family weight ≥ 5).

## How the sample was reached
71 texts in MGH Auctores antiquissimi 1–15 recorded blind; 46 with ≥ 2 sigla and a top (36 explicit, 5 implied, 5
undecided). Source: vision-model reading of the dMGH page images (A8), parser v4 (A10). Located: 42. Screen (≥ 4
analysis witnesses cited ≥ 30 times, ≥ 200 units with a shared departure): 28. Proof page (precision ≥ 0.85, recall
≥ 0.60, blind gold from the image): passed by 21. Of these, 14 have an explicit top and are not pilots: the primary set.
Repair round (A10): 9 texts failing the proof page got one parser repair, diagnosed on 18 development pages only, and
fresh proof pages. 6 passed, 3 did not (Acta synhodorum, Aldhelm De virginitate prosa, Symmachus epistulae).
Excluded at the proof page after the repair: also Fortunatus carmina (pilot, failed the re-score) and three Aldhelm
texts newly passing the screen.

## Primary result (explicit tops, 14 texts; all k_E = 2 except Romana, k_E = 1)
| | n | 95% CI (Clopper–Pearson) |
|---|---|---|
| **supported** | 4 | 0.08–0.58 |
| **contradicted** | 4 | 0.08–0.58 |
| **partly silent** | 6 | 0.18–0.71 |
No k_E ≥ 3 group, so no test (as pre-registered).

| Text (editor) | Editor's top | Verdict | Evidence (stable families, weight) |
|---|---|---|---|
| Sidonius, carmina (Luetjohann, filiation by Leo) | M \| TCFP | supported | TCFP 32 |
| Eugippius (Sauppe) | L \| VMD | supported | VM 214, VMD 51 |
| Marcellinus Comes (Mommsen) | S \| TUR | supported | UR 86, TUR 43 (0.91) |
| Victor Vitensis (Halm) | BLV \| R | supported | BV 264, BLV 149, BLRV 58 |
| Ennodius (Vogel) | B \| VLTCOP | **contradicted** | BVL 211 and BVLT 71 join B to V L (T) |
| Vita Martini (Leo, Fortunatus) | PTF \| NBLMS | **contradicted** | PN 75 crosses both branches |
| Paulus, Historia Romana (Droysen) | HVBM \| NZ | **contradicted** | NZM 48 joins M to N Z |
| Prosper (Mommsen; admitted under A8c) | AORV \| MY | **contradicted** | AC 108, ACV 21 cross AORV |
| Claudian, carmina maiora (Birt) | VPGM \| CBAJR | partly silent | VP 362, BA 191, nothing at the top |
| Claudian, Raptus (Birt) | FVS \| CABDW | partly silent | CB 22 only |
| Historia Brittonum (Mommsen) | HKPQDCGL \| ZMN | partly silent | MN 686, PQDCGL 286, HK 36 |
| Isidore, Chronica maiora (Mommsen) | 12 \| 7 witnesses | partly silent | low-level pairs only |
| Isidore, Historia Gothorum (Mommsen) | P \| BNHOXYRMA | partly silent | MA 80, XY 13 |
| Jordanes, Romana (Mommsen) | HPVL | partly silent | HPV 136, HP 74 |

## Sensitivity (as pre-registered and per A8c/A10)
- Without the four repaired texts (Isidore CM, Eugippius, Victor, Ennodius): 10 texts, 2 supported / 3 contradicted / 5 silent.
- Without Prosper: 13 texts, 4 / 3 / 6.
- Prosper, Paul and the empty proof pages (Eutropius, Paul, Ausonius: nothing printed, nothing parsed) say little about extraction there.

## Reported separately
Pilots (used to build the method): Sidonius epistulae **supported** (MTCFP 365); Getica **contradicted** (Mommsen's
I+II branch HPVLAOB crossed by AOBXYZ 7, 0.93: A with O B and the third order; consistent with Grillone 2017);
Variae **contradicted** (KEF 197 crosses both classes); Eutropius **contradicted** (LOD 123 crosses Droysen's GFD|LO).
Implied top (A3): Eugenius Dracontiana contradicted (MFP 10). Undecided (A3): Ausonius contradicted (BZ 29), Avitus
contradicted (NL 121).

## Reading
Of fourteen explicit top divisions, the editor's own apparatus confirms four, contradicts four with a stable shared
departure that crosses the division, and is silent on six: the lower families are stable but nothing groups the top.
The contradictions are of one kind: a family of manuscripts that the editor places on opposite sides of his top split
shares many departures (Ennodius B with V L: 211 units). They can mean contamination or conjunctive errors read as
kinship; the method does not decide which, but it names the units. Silence is the commonest outcome: the evidence for
the top of these stemmata is not in the apparatus, which is Bédier's suspicion in its weak form. Bipartition as such
is not tested: every explicit top here is two-way.

## EXPLORATORY (not pre-registered, run after the classification): editors' tops against random tops
code/historic/bedier_null.py; data/historic/bedier/null_exploratory_primary.json. For each text the stable families of
the frozen run are kept, and every pair of disjoint branches with the editor's branch sizes (drawn from the analysis
witnesses) is scored by the frozen rule; enumeration is exact. 13 texts (Romana has one branch only).
| | editors | expected under random tops | one-sided p |
|---|---|---|---|
| supported | 4 | 0.94 | 0.005 |
| not contradicted (supported or silent) | 9 | 2.27 | 8 × 10⁻⁷ |
| score (supported 2, silent 1, contradicted 0) | 13 | 3.21 | 1.4 × 10⁻⁵ |
Random tops are contradicted 40–100% of the time (median ≈ 0.86); editors' tops are contradicted in 4 of 13. So the
editors' divisions carry real signal: far better than chance, even where the apparatus does not prove them.
Caveats: rooting uses the editor's own text, which he built with his stemma in mind, so the comparison is not free of
his judgement; texts are treated as independent; exploratory, decided after the classification.

## EXPLORATORY: why two-way? Two checks (code/historic/bedier_why.py; data/historic/bedier/why_exploratory.json)
**Merge check.** "Joins" = shared departures linking ≥ 2 of the data's top blocks (maximal stable families + loose
witnesses). For the editor's split: share of joins inside a branch vs across the split, against every (or 20,000
random, for > 50,000) splits with the same branch sizes.
- Supported texts (data top = editor's top, k_D = 2): the branches are single blocks, so every join is across; these
  weak crossing groups (10–591 units) are the noise floor that even a well-supported division carries.
- Editor's split captures joins far above random: Isidore Chron. mai. 0.44 vs 0.28 (p < 0.001), Variae 0.54 vs 0.13
  (p < 0.001), Prosper 0.38 vs 0.12 (p 0.018), Getica 0.49 vs 0.26 (p 0.03).
- Not above random: Vita Martini, Paul, Raptus, Hist. Britt., Eutropius (p 0.10–0.36); Claudian maior below random
  (0.17 vs 0.37).
- Even where it is real, the signal for the editor's top is at most half of all joins: as many shared departures cross
  the split as support it. The editor picks one dichotomy from alternatives of comparable weight.
**Choice check** (both branches ≥ 2 witnesses, 10 texts): units where one branch departs as a whole and the other is
intact; the editor's text then follows the intact branch. It leans on one branch in every text: Vita Martini 90%
(FPT), Paul 83%, Raptus 83%, Isidore CM 80%, Claudian 78%, Eutropius 78%, Getica 71%, Variae 71%, Hist. Britt. 66%,
Prosper 58%; not always the larger branch (Vita Martini, Variae, Eutropius follow the smaller one).
**Reading.** Consistent with both of Bédier's mechanisms: the evidence at the top is weak and balanced, and the two-way
top leaves the editor free to follow a preferred branch, which he does about four times in five. Not a proof of motive;
exploratory, decided after the test.
