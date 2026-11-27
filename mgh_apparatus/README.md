# MGH *Auctores antiquissimi*: critical apparatus in tables

The critical apparatus of 42 texts edited in the *Monumenta Germaniae Historica, Auctores antiquissimi*
(volumes 1–15, 1877–1919), read by machine from the page images and set out as plain tables: 122,054
apparatus entries, 150,947 readings, each with the manuscripts the editor gives for it, its page image and the
transcribed line of the apparatus it comes from. Made for the series test in *Testing stemmata against the
evidence* (§4.3), and released for anyone who wants to query these apparatus without re-reading the volumes.

Metamorphoses of Civilization Collection, 2026. Tables: CC BY 4.0. The MGH volumes are in the public domain;
page images: dMGH, https://www.dmgh.de/.

## Files

| File | Content |
|---|---|
| `apparatus/<text_id>.csv` | one row per reading: `entry` (apparatus entry number within the text), `volume` (dMGH volume id), `image_leaf` (dMGH page image), `text_line` (line of the edited text the entry refers to, as printed), `lemma` (when printed before `]`), `reading`, `witnesses` (sigla carrying the reading), `apparatus_line` (the transcribed apparatus line) |
| `texts.csv` | one row per text: work, editor, year, volume, first and last page image, number of entries and readings, witnesses cited (by frequency), the editor's top division as stated in his preface |
| `sigla.csv` | the sigla of each text as the editor's preface lists them, with manuscript, date and role |
| `proof_pages.csv` | every proof page used to measure the extraction: a page drawn at random, its apparatus transcribed from the image, scored against the tables (groups on the page, recovered, extra, precision, recall), by extraction round |
| `../data/bedier/ocr_llm/` (`page_transcriptions/` in the Zenodo zip) | the 4,091 page transcriptions the tables were extracted from, one file per page image |

Rebuild: `python3 code/historic/export_apparatus_tables.py <data dir> mgh_apparatus`.

## How the tables were made, and how far to trust them

The page images (dMGH) were transcribed by a vision model (Gemini, `google/gemini-4-argon`, temperature 0;
prompt in `code/historic/bedier_llm_read.py`). A single extractor for MGH-style negative apparatus
(`code/historic/mgh_rooted_extract.py`) then split each apparatus line into entries and readings and kept the
sigla. This is automatic extraction, not a proofread digital edition:

- `witnesses` keeps the sigla of the edition's main witnesses. Sigla of excerpts, later hands and witnesses
  outside the editor's stemma are often dropped (*CAL* can appear as *CA*); the full printed form is in
  `apparatus_line`.
- Readings are sometimes cut at the wrong place, or merged with the next entry, where the apparatus is dense.
  Check the `apparatus_line`, and the page image for anything that matters.
- Precision and recall of the shared groups were measured on proof pages (`proof_pages.csv`). The texts used in
  the series test passed a page drawn at random with precision at least 0.85 and recall at least 0.60.
- A blind hand check of 60 readings on the page images (40 behind the paper's four contradicted top divisions,
  20 from other texts) found the witnesses correctly recorded in all 60; one reading was misspelled.

Cite as: Metamorphoses of Civilization Collection (2026), *stemma-mwis: code and data for "Testing stemmata
against the evidence"*, https://doi.org/10.5281/zenodo.24405388, folder `mgh_apparatus/`.
