# Transcription spec: Licciardello 2018, Passio Donati III (BHL 2294), apparatus pp. 264–313

COPYRIGHT: study copy. Output stays in a local, non-public folder (never committed, never shared).

Page images: a local folder of the edition's pages, one image per page (pNNN.jpg)
Each page = edited text (top, with chapter Roman numerals and sentence numbers) + apparatus (bottom, small type).
Read the images yourself with the Read tool (zoom by reading again if needed; you may crop with
`sips`/python PIL into your own scratch dir to enlarge the apparatus band). Do NOT use OCR.

## Output: one JSON file per page
`<OUTDIR>/pNNN.json`:
```json
{
 "page": 270,
 "chapter_at_top": "IV",            // chapter in force at the top of this page's text
 "coverage": [                       // every bracketed absence note, e.g. "[om. Nn, U1]" before sentence 6
   {"chapter": "IV", "sentence": 6, "wits": ["Nn","U1"], "raw": "[om. Nn, U1]"}
 ],
 "units": [
   {"chapter": "IV", "sentence": 6, "lemma": "humani generis",
    "readings": [
      {"text": "generis humani", "wits": ["L2","LAZZ"], "type": "order"}
    ],
    "raw": "humani generis: generis humani L2, LAZZ",
    "uncertain": false}
 ],
 "other_notes": []                   // anything else about witness presence (inc./des./mutilo/deest/lacuna/fragment), verbatim with chapter.sentence
}
```

## Rules
- One unit per lemma (each "lemma: variant W; variant W ~" block; `~` separates units, a bold/number starts a new sentence).
- Witness tokens exactly as printed, normalised: Latin sigla as `F14`, `L1`, `U1`, `SO`, `Gb`, `Nn`; editions as `BUON`, `DRAG`, `LAZZ`, `MOMB`.
- Group sigla in Greek, written as these ASCII names: kappa, lambda, mu, nu, xi, omicron, pi, rho, sigma, tau, upsilon, phi.
  (Careful: in the apparatus `o` = omicron group, `v`/`ν` = nu, `υ` = upsilon, `ϱ`/`ρ` = rho, `x`/`κ` = kappa.)
  With exceptions: write `"nu(-rho,-V3)"` for "ν (tranne ϱ, V3)"; `"kappa(-H)"` for "κ (tranne H)".
- Corrections: "vocantes corr. in invocantes Li" → Li read *vocantes* first, then *invocantes*.
  Encode the hands as `Li^ac` (first hand) and `Li^pc` (corrected): put `Li^ac` under the reading it had first and
  `Li^pc` under the reading after correction (if that is the lemma/critical text, add a reading with `"text": "<lemma>"`).
  Same for "agg. in interl." (addition above the line: `R3^ac` lacks it, `R3^pc` has it), "segue una parola erasa", "in marg." etc.
- `type`: one of `spelling` (orthography only: ae/ę/e, ti/ci, h-, y/i, doubled/single consonants, assimilation,
  aruspici/auruspici-type, inflection-neutral spelling of names), `word` (different word/form), `order`,
  `omission` (om.), `addition` (segue X / add.), `other`. When unsure between spelling and word, use `word` and set `"uncertain": true`.
- `uncertain: true` whenever a siglum, a digit (1/I, 3/5/8, 4), or a group letter is not certain on the image. Never guess silently.
- Skip nothing: every entry on every page, including editions-only entries.
- Include the critical-text chapter numbering exactly; if a page starts mid-chapter, look at the previous page to know the chapter.
- When done, write `<OUTDIR>/DONE_<first>-<last>.txt` with: pages done, number of units, number marked uncertain.
- Report back only counts and any systematic difficulties. Do not paste apparatus text into your reply.
