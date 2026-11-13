# Brief: transcribe the critical apparatus of one printed page (gold standard, blind)
You read page images of 19th-century MGH editions (Latin, negative apparatus at the foot of the page). For EACH page in
your batch file, transcribe every apparatus entry that cites at least one manuscript siglum, as printed:
  {"line": <line number printed in the apparatus, or null>, "lemma": "...", "readings": [{"text": "...", "sigla": "XY2Z"}]}
- Write sigla exactly as printed, in order, including hands/superscripts (write superscript 1 as "1", 2 as "2", e.g. "P1",
  "C2"), Greek letters as Greek, and abbreviations like "Bern." as printed. Separate sigla by nothing ("LMTCFP") unless
  printed apart; if a siglum has more than one letter, put spaces between sigla ("A Bern. R").
- If the lemma itself carries sigla (lemma] SIGLA), put them in "lemma_sigla".
- Include every entry on the page, also single-siglum ones; skip notes that cite no siglum, and skip the edited text.
- Use the Read tool on the image (it displays it); zoom by cropping is not available, so read carefully, and if a siglum is
  not legible write "?" for it. Do NOT open any other file of the project (no OCR text, no units.json, no results): this
  must be a reading of the image alone.
Output: one JSON file per batch, at the path given, {"pages":[{"text_id":..,"page":..,"entries":[...],"notes":"..."}]};
validate with python3 -c "import json;json.load(open(PATH))". Report one line per page: text_id | n entries | legibility.
Each page also has "bands_top_to_bottom": overlapping horizontal crops of the lower page at full resolution (y from 1450 to
~3300 px, 620 px high, 170 px overlap). Read the apparatus from the bands (use the full image only to see where the
apparatus starts); beware of entries repeated in the overlap of two bands: transcribe each printed entry once.
