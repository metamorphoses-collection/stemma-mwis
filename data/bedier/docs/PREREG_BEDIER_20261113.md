---
type: pre-registration
project: priscus_excerpta (paper: Reconstructing manuscript genealogies…, new Section "Bédier's question")
date: 2026-11-13
status: FROZEN when committed; any later change goes in a dated "Amendments" section below, never by editing above it
code: code/historic/mgh_rooted_extract.py (v2, unchanged), code/historic/rooted_mwis.py, code/historic/bedier_top.py (to be written; scoring only)
---

# Bédier's question on the editors' own evidence: pre-registration

## Question
In 1928 Bédier observed that almost all published stemmata divide into exactly two branches at the top, and
suspected that this reflected the editors, not the traditions. Here is the question for each stemma: **does the
apparatus the editor printed support the editor's own top division?** We count shared departures from the edited
text, as Lachmann's rule requires, with exact rooted MWIS and bootstrap. Our second question: are the
two-branch tops supported any differently from tops with three or more branches?

This measures consistency of each stemma with the evidence its editor printed, not the true history.
Either answer is a result. "The fathers' bipartitions hold up on their own evidence" and "they do not" are both
reportable; the protocol below is fixed so that neither can be chosen after the fact.

## Sampling frame (closed)
Every text edited in *Monumenta Germaniae Historica, Auctores antiquissimi*, vols 1–15 (1877–1919), using the
Internet Archive scans (University of Toronto copies: monumentagermani01geseuoft … monumentagerman15geseuoft;
the gese/Google copy where no uoft copy exists). Proof-checks use dMGH page images.

A **text** is a work with its own manuscript tradition, apparatus and stated filiation. A volume holding several works
contributes several texts. A work edited in several sections with separate witness sets (e.g. the Variae) is one text,
unless the editor gives separate stemmata.

## Inclusion, applied in this order and logged with the reason for every exclusion
1. The editor states a filiation or a division into classes or families of manuscripts in the preface (explicit
   diagram or explicit verbal statement). Otherwise the text is excluded as "no stemma".
2. The apparatus is a negative or positive apparatus citing manuscripts by sigla. Otherwise: "no sigla apparatus".
3. At least 4 witnesses of the stemma are cited in at least 30 extracted departures each (the **analysis witnesses**).
   Otherwise: "too few witnesses".
4. At least 200 extracted units with a shared departure (≥ 2 analysis witnesses). Otherwise: "too thin".
5. Extraction quality on one proof page (below): precision ≥ 85% and recall ≥ 60%. Otherwise: "OCR".

## Recording the editor's top division (blind)
For each included text, before its analysis is run, the editor's top division is written to
`data/historic/bedier/editor_tops.json`, with page references to the preface. A top division is the set of branches
directly below the archetype, restricted to analysis witnesses. A branch with no analysis witness is dropped; a
single witness can be a branch. k_E is the number of branches. Rules:
- If the editor gives several archetypes or allows contamination at the top, record the main tree and flag it.
- If the editor is undecided between options, record the text as "undecided" and report it separately.
- Indirect witnesses (excerpts, testimonia, early editions) are not analysis witnesses.
The file is committed before the analysis of those texts is run; the commit hash is cited in the results.

**Pilot set (not blind).** Getica, Eutropius, Sidonius (letters), Variae, Fortunatus were analysed on 12 Nov before this
protocol. They are scored by the same rules but reported separately and excluded from the primary counts.

## Extraction and analysis (unchanged code)
- Extraction: mgh_rooted_extract.py v2 on the OCR of the text's page range. Sigla, ranges and skip regions are set from
  the preface and page heads, before any results are seen. Per-text settings go to `data/historic/bedier/settings.json`.
- Groups equal to all analysis witnesses are dropped (editorial conjectures against all witnesses). Non-analysis
  witnesses are removed from every group (`--drop`).
- Rooted MWIS (rooted_mwis.py), minimum weight 2, B = 200 bootstrap replicates over units, seed 1.
- A **stable family** is a family in the exact optimum with bootstrap stability ≥ 0.90.

## Proof page (one per text, chosen before extraction)
The page is drawn by Python `random.Random(<volume>*1000 + <text index>)` uniformly among the text's apparatus pages,
and recorded in settings.json before extraction. On the dMGH image of that page, every printed group of ≥ 2
first-hand sigla sharing a reading is counted; parsed groups on the same page are matched. Precision and recall are
computed as in data/historic/mgh/proof (gold_groups.json).

## Primary outcome: classification of the editor's top division
Each top branch B_i with ≥ 2 analysis witnesses is classified:
- **supported**: B_i is a stable family;
- **contradicted**: some stable family of weight ≥ 5 crosses B_i (overlaps without nesting);
- **silent**: neither.
Single-witness branches are not classified.

The text's top division is **supported** if every classified branch is supported; **contradicted** if any branch is
contradicted; otherwise **partly silent**.

Reported: counts of supported / partly silent / contradicted, separately for k_E = 2 and k_E ≥ 3. These are reported
with exact (Clopper–Pearson) 95% intervals, and with no significance test unless both groups have ≥ 8 texts. In that
case the test is Fisher's exact test on supported vs not, two-sided.

## Secondary outcomes (descriptive)
- k_D: the number of top-level blocks in the data. These are the maximal stable families, plus analysis witnesses in no
  stable family. Also the bootstrap distribution of k_D, using maximal families of each replicate's optimum.
  This is reported with the caveat that missing resolution raises k_D, so k_D > k_E alone is not evidence against
  the editor.
- For each contradicted text: the crossing family, its weight, stability and five sample loci. These are checked on page
  images where digitised, as was done for the Ambrosianus and S. Marco 554.
- Known bias, stated in advance: rooting by the edited text makes the branch the editor follows look silent. Silence is
  therefore reported as silence, never as contradiction.

## What would change the paper
Nothing is decided in advance about the wording. The section reports the table of all included texts, the exclusion
log, the pilot set separately, and the classification counts.

## Amendments
### A1 (2026-11-13, before any confirmatory text was extracted or any editor's top recorded): OCR source
CJ asked for the volumes to be OCR'd "with utmost precision", with several techniques if needed. The OCR source is
therefore chosen by a fixed benchmark rather than by assumption.
- Candidates: (a) Internet Archive OCR (djvu.txt); (b) Apple Vision on dMGH 300 dpi images (accurate level, language
  correction off; Vision has no Latin model, so three language settings are benchmarked: it-IT, en-US, default); (c) Tesseract `lat` on the same images; (d) combinations of these at the level of parsed
  groups: union, and groups found by at least two engines.
- Benchmark: the existing 200 gold groups on 12 pages (data/historic/mgh/proof/gold_groups.json), scored by
  proof_score.py. These pages come from the pilot editions, so the benchmark sees no confirmatory text.
- Rule: choose the candidate with the highest recall among those with precision ≥ 90%; break ties by precision. The
  scores and choice are recorded here before confirmatory extraction starts.
- All page-level OCR outputs are kept per leaf (`@@LEAF n` markers), so every extracted unit traces to a page image.

### A2 (2026-11-13, before any confirmatory text was extracted or any editor's top recorded): OCR benchmark result
Scores on the 200 gold groups (12 pages), with the v3 parser. v3 adds one OCR repair: a capital I after a siglum,
where I is not a siglum, is read as the hand mark 1. The Internet Archive score is unchanged by it.
| Candidate (as listed in A1) | Recall | Precision |
|---|---|---|
| (a) Internet Archive OCR | 169/200 = 84.5% | 90.4% |
| (b) Apple Vision, whole page (it-IT = en-US = default, byte-identical) | 61% | 94% |
| (b') Vision, 2x upsample / strips 240 px 2x (edge rule) / strips 360 px 3x | 58% / 65% / 60% | 96% / 94% / 96% |
| (c) Tesseract lat (tessdata_best) psm 6 / psm 4 | 73% / 71% | 91% / 92% |
| (d) union of IA, Vision strips, Tesseract | 94% | 84% |
| (d) groups found by >= 2 of {IA, Vision strips, Tesseract} | 82% | 98% |
Vision reads very accurately but drops short, isolated apparatus entries, and its line order mixes entries printed side by
side. Combinations designed after seeing these scores are not eligible.
**Primary OCR source (rule of A1): (a) Internet Archive OCR.**
**Sensitivity analysis, added now:** every confirmatory text is also extracted from the >= 2-engine vote (precision 98%).
The classification is reported under both. Any text whose classification differs between the two is flagged and listed.
For this, the Studio OCRs all volumes with Vision strips (both geometries) and Tesseract lat psm 6 and 4. All engines'
per-leaf outputs are kept (data/historic/bedier/ocr_apple/), also as a corpus in their own right.
Scripts: code/historic/ocr_bench/bench.py, combo2.py.

### A3 (2026-11-13, before any confirmatory text was extracted): how recorded tops are classified
The blind step produced readings of four kinds. They are fixed now, before extraction:
- **explicit**: stated in the editor's words, or in a printed stemma read on the page image. These enter the primary counts.
- **implied**: the agent's grouping follows from the editor's remarks (e.g. "H ex G, B ex H", "A praecedenti persimilis"),
  but the editor states no top division. These are reported separately with the undecided texts, never in the primary counts.
- **undecided**: the editor leaves the top open or offers it tentatively. Reported separately (as in the main protocol).
- **none**: no filiation stated. Excluded as "no stemma".
Where an agent's reading and a printed stemma disagree, the stemma read on the image is recorded (Vita Martini, Eugenius
carmina). Every decision is in data/historic/bedier/tops/VERIFICATION.md. Implied: Vita Hilarii, Vita Germani, Expositio
fidei, Aldhelm Carmina ecclesiastica, Eugenius Dracontiana (rests on a garbled footnote). Undecided: Chronographus 354,
Victorius, Polemius, Ausonius, Avitus carmina.

### A4 (2026-11-13, before any confirmatory extraction): extraction conventions
Fixed mechanically from the frozen records (tops/*.json), before any text is extracted:
- **Sigla form.** The extractor (v3) reads sigla that are single capital Latin letters. Witnesses with other sigla (Greek,
  lower-case, several letters, or a letter + digit naming a separate manuscript) cannot be read reliably from OCR. They are
  treated as non-analysis witnesses. A text left with fewer than 4 analysis witnesses is excluded as "too few witnesses
  (siglum form)".
- **Digits.** Where the edition uses letter + digit for distinct manuscripts (Aldhelm, Vita Radegundis: P1–P4, S1, S2,
  X1, X2, Y1, Y2 …), extraction runs with `--digits drop`: such tokens are dropped, never read as hands. Elsewhere digits
  are hands (X, X1 = first hand kept; X2+ = corrector, dropped), as in the pilot.
- **OCR character fixes:** '0:O' if O is a siglum; '8:S' if S is a siglum; nothing else.
- **Ranges.** The text_lines of the record. Where one text is printed inside another (Aldhelm's Aenigmata inside De metris;
  Eugenius' Dracontiana interleaved with Dracontius; Eutropius and Paulus on the same pages), each is extracted from its
  own range with only its own sigla. Overlaps are logged; the scorer's sigla filter does the rest.
- **Proof page:** the printed page is drawn as in the protocol, with volume code n for AA n, and 61 or 62 for AA 6,1 or 6,2.
  The text index is the record's position in its tops file. The draw is among page numbers found in the text's OCR range.
- **Witnesses passed to the scorer:** the record's direct sigla of single-capital form, in record order. Top branches are
  the record's, restricted to those.

### A5 (2026-11-13): locating the proof page in the OCR (measurement, disclosed in full)
Rule 5 needs the OCR lines of each drawn page. Two page-number heuristics proved unreliable in the Internet Archive OCR
(first: 2/21 pass, with 7 pages not found; second: results moved in both directions). They were abandoned, not tuned.
The rule adopted is content alignment (code/historic/bedier_align.py):
- Vision reads the first lines of edited text on the page image and on the next page.
- These are fuzzy-matched in the OCR. The window runs from this page's match to the next page's match, 30–250 lines later.
- A page is located only if both matches score ≥ 0.75; an unlocated page fails rule 5.
A bug was also found and fixed: units carried the line where their apparatus block starts, not their own line. The
extractor now also records `entry_line`. All 46 re-extractions are identical in every lemma, group and reading.
**Result with the primary source (Internet Archive OCR):** 2 of 21 texts pass rule 5 (Eugippius; Fortunatus carmina,
pilot). 6 pages are not located (Paulus' page is a Paeanius page; Eutropius' page is Greek). 13 fail on recall (precision
is mostly ≥ 0.8). The cause is the OCR of sigla, e.g. Romana p. 19: HPL → "IIPI.", HPV → "IJPV".
Scores: data/historic/bedier/proof_scores.json; gold transcriptions: data/historic/bedier/proof_gold/.

### A6 (2026-11-13, approved by CJ; before any classification was run): source after the rule-5 failure
The primary source failed rule 5 on 19 of 21 texts (A5). As CJ approved, the source becomes the per-leaf multi-engine
OCR (Studio, dMGH 300 dpi). It was designated in A2 as the sensitivity source. The gold transcriptions, the proof pages,
rules 1–5 and the scorer stay exactly as they are.
- **Combination rule**, chosen by the A1 rule (highest recall with precision ≥ 90%) on the 12 *benchmark* pages, not on
  the proof pages:
  | per-leaf candidate | recall | precision |
  |---|---|---|
  | Vision strip2 ∪ strip3 | 72% | 94% |
  | Tesseract psm6 ∪ psm4 | 74% | 91% |
  | Vision ∩ Tesseract | 58% | 97% |
  | Vision ∪ Tesseract | 88% | 89.8% (fails) |
  | **≥ 2 of the 4 variants agree** | **81%** | **92.6%** |
  | ≥ 3 of 4 | 58% | 97% |
  Selected: a group (any size, single witnesses included) is kept on a leaf with the multiplicity at which at least two
  of the four variants (strip2, strip3, tess, tess4) report it.
- **Units.** Each variant is extracted per leaf with the text's settings (sigla, digits, OCR fixes). The voted groups of
  a leaf become units of one group each. Family weights equal those of the voted counts. The bootstrap resamples these
  group-units.
- **Leaves of a text.** The text's printed page range (first to last page found in its IA range), mapped to leaves by the
  volume offset measured on the proof page. The mapping is checked by Vision's page number on the first and last leaf.
- **Inclusion (rules 3–4) is recomputed** on the new source, and rule 5 is scored per leaf (no alignment needed).
- The IA-based extraction stays as the record of A5. It is not used for classification.

### A7 (2026-11-13, before any proof page or classification on the A6 source): which leaves belong to a text
A6 mapped a text to leaves through its printed page range. A first run showed that the page-number runs found in the IA
OCR cover only part of long texts (Ennodius: pages 51–67 of several hundred; Getica: 91–137), and that the volume part
(e.g. AA 5,2 for Ausonius, 4,2 for Krusch's Vitae) was not identified. The leaves are therefore located by content,
with the A5 method (code/historic/bedier_a6.py):
- Four lines of edited text (≥ 6 words, ≥ 80% lower-case words, no sigla) are taken after the start of the text's IA
  range, and four before its end.
- They are fuzzy-matched against the strip2 OCR of every leaf of every part of the volume; the end is searched after the
  start, in the same part. Both matches must reach 0.75, or the text is not located and is excluded ("not located").
- A volume counts as OCRed only when its phase-2 output is present (≥ 50 leaves).
First run on AA 2–7: 15 texts located; Expositio fidei (implied) and Avitus carmina (undecided) not located.

### A8 (2026-11-13, approved by CJ: "strengthen ocr"; before any classification, before any text was read with the new source): vision-LLM page reading
**Why.** On the A6 source, proof pages of AA 8–12 were scored before this amendment: Claudian passes; Sidonius (letters,
poems), Cassiodorus' Variae, Isidore (Historia Gothorum, Chronica maiora) fail; the Roman synods' proof leaf fell outside
the located leaves. The failures are apparatus reading (long sigla runs lost, short entries dropped), not the method.
These scores motivated the amendment; nothing below was chosen or tuned on any proof page.
- **Source.** Each leaf of a located text is transcribed in full from the dMGH 300 dpi image by a vision model
  (google/gemini-4-argon via OpenRouter, temperature 0; prompt frozen in code/historic/bedier_llm_read.py), then parsed by
  the unchanged extractor with the text's frozen settings (sigla, digits, OCR fixes). Single source, no vote.
- **Selection, on the 12 benchmark pages only** (gold of A2, same scorer):
  | source | recall | precision |
  |---|---|---|
  | A6: ≥ 2 of 4 engines | 81% | 92.6% |
  | **vision-LLM page reading** | **89%** | **97%** |
- **Location (A7, two mechanical additions, blind to outcomes).** (i) If no anchor line of edited text is found within
  600 IA lines of a range end (tabular endings, e.g. Isidore's Chronica maiora), the scan is widened to the whole range.
  (ii) After the located end leaf, following leaves are added while their running head matches the end leaf's
  (ratio ≥ 0.8, at most 15 leaves): texts ending in lists without prose (the synods' subscriptions) were cut short.
  Neither can move a text located before. Not located, hence excluded: the AA 9 chronicles (garbled IA OCR; the dMGH
  table of contents was checked and its titles are offset, so it is not used), Avitus carmina, Expositio fidei.
- **Proof pages.** The pages drawn under A4 and their blind gold stay; they are scored once on the new source. Gold was
  read by Claude agents, the source by a Gemini model (different families, independent errors). Texts newly passing the
  screen get a proof page drawn by the A4 rule and transcribed blind before scoring.
- **Inclusion (rules 3–4) is recomputed** on the new source. Classification runs once, after all proof pages are scored.
- **Stored.** Transcriptions (text only, no images) in data/historic/bedier/ocr_llm/<vol>/, with USAGE.tsv (tokens,
  cost per leaf). The A6 vote stays as the record and as a sensitivity source.
- **A8b (same day, before these texts were read).** If a text's start is located (≥ 0.75) but its end is not, the text
  ends at the leaf before the next section boundary of the dMGH table of contents (page labels only, titles unused,
  cached in data/historic/bedier/dmgh_toc/). Recovers Polemius (undecided), Chronica Gallica (explicit) and Avitus
  carmina (undecided). Still not located: Chronographus, Liber generationis, Victorius, Expositio fidei, and Prosper
  (start match 0.749, below the threshold, left unchanged).
- **A8c (same day; deviation, decided by CJ: "i want prosper"; before Prosper was read or screened).** Prosper's start
  matched at 0.749 (leaf 357 = p. 346, Prosper's preface; checked on the running head), just below 0.75. It is admitted
  with that start and the A8b end (next dMGH section boundary: leaves 357–512). It must pass the screen and a proof page
  drawn by the A4 rule like any text. All results are reported with and without Prosper.

### A9 (2026-11-13, before any proof page of these texts was read or scored): proof pages of texts newly passing the screen
On the A8 source the screen (rules 3–4, frozen scorer) includes three explicit texts that had no scored proof page:
Claudian's Raptus Proserpinae, Marcellinus Comes, Prosper (A8c), and the undecided Ausonius. Their proof pages were drawn
in settings.json before extraction (A4 rule), as printed pages located through the IA page chain.
- Raptus (p. 379 → leaf 613) and Ausonius (p. 221 → leaf 286) fall inside the located leaves (dMGH page labels): kept.
- Marcellinus (p. 388) and Prosper (p. 53) fall outside their texts (the IA page chain was wrong there). They are
  redrawn with the same seed, random.Random(volume*1000 + index), uniformly among the text's located leaves that carry
  an apparatus (a transcribed line containing "]"): Marcellinus leaf 67 (seed 11000), Prosper leaf 427 (seed 9005).
- Gold for the four pages is transcribed blind from the dMGH image (brief of proof_gold/brief.md, reader sees no OCR,
  no extraction) and scored once with bedier_proof.py.

### A10 (2026-11-13, approved by CJ: "Yes, one repair round"; before any diagnosis): one repair round for texts failing rule 5
On the A8 source, eight included explicit texts fail their proof page (Variae, Acta synhodorum, Aldhelm De virginitate
prosa, Isidore Chronica maiora, Eugippius, Victor Vitensis, Symmachus epistulae, Ennodius), and the undecided Avitus
carmina (added at CJ's request, "and fix the pilots + ausonius and avitus"). One repair round, then no more:
- **Pages drawn now, before any gold** (proof_gold/a10_draws.json): for each of these nine texts,
  random.Random("A10-" + text_id).sample(apparatus leaves of the text minus the old proof leaf, 3): the first two are
  **development pages**, the third is the **fresh proof page**. The old proof pages are retired: their gold and parsed
  groups are not inspected and not used for development.
- **Gold**: all development and fresh proof pages are transcribed blind from the image (brief of proof_gold/brief.md).
  The fresh-proof gold is stored separately and not opened before the parser is frozen.
- **Diagnosis and fix**: on the development pages and the 12 benchmark pages only. Fixes go into the parser or the
  per-text frozen settings (sigla, OCR fixes), apply to every text, and are listed in this file before scoring.
- **Then**: every text is re-extracted with the frozen parser; the fresh proof pages of the nine are scored once; the
  proof pages of all other texts (passed under A8/A9) are re-scored once with the new parser, and a text that now fails
  is excluded. The screen is recomputed, and the classification runs once.
- Disclosed in the results: the A8 failures, this round, and the classification with and without the repaired texts.
- **A10 fixes (parser v4, frozen before any fresh proof page or re-score was looked at).** Diagnosed on the 18
  development pages (gold_dev_1/2.json) only. All apply to every text:
  1. "lemma SIGLA:" — a sigla token ending in ':' before any group is the lemma's own sigla, not a departure (Victor, Eugippius).
  2. Capitals outside the text's sigla inside a sigla run of ≥ 2 letters are dropped, the rest kept ("is BZNEFA" → BEFNZ
     in the Variae, where A is not an analysis witness); before, the whole run was rejected.
  3. Lower-case letters: settings 'lower' lists the single lower-case sigla of the text's blind record (Victor p, Avitus
     g n p, …). After a capital such a letter is a separate non-analysis siglum (dropped, capital kept), any other
     lower-case letter is a superscript (that siglum dropped, as the frozen scorer did); standalone it ends an entry.
     The proof scorer applies the same rule (gold_groups 'lower').
  4. Footnote keys "a) b) …" start entries and mark apparatus lines (Aldhelm).
  5. Interlinear apparatus without lemma or line number: lines with ≥ 2 "word SIGLA," pairs, lines starting
     "om. SIGLA", and lines "word SIGLA" (Acta synhodorum).
  6. A group repeating a siglum already cited in the current entry starts a new entry (was: whole entry dropped).
  Development pages: recall 0.72 → 0.92 (batch 1), 0.86 → 0.955 (batch 2); precision 0.93 → 0.97, 0.84 → 0.90.
  The 12 benchmark pages (vision reading): recall 89% → 94%, precision 97% → 95%.
