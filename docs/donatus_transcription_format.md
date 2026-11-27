# Donatus apparatus: transcription format

The Donatus matrix is built by `code/donatus/build_matrix.py` from a transcription of the negative apparatus of
P. Licciardello's edition (*La «Passio» di san Donato vescovo di Arezzo*, SISMEL, 2018, pp. 264–313), which is
not redistributed. The transcription is one JSON file per printed page, `pNNN.json`:

```json
{
 "page": 270,
 "chapter_at_top": "IV",
 "coverage": [
   {"chapter": "IV", "sentence": 6, "wits": ["Nn", "U1"], "raw": "[om. Nn, U1]"}
 ],
 "units": [
   {"chapter": "IV", "sentence": 6, "lemma": "humani generis",
    "readings": [{"text": "generis humani", "wits": ["L2", "LAZZ"], "type": "order"}],
    "raw": "humani generis: generis humani L2, LAZZ",
    "uncertain": false}
 ],
 "other_notes": []
}
```

- `chapter_at_top`: the chapter in force at the top of the page. `coverage`: every bracketed absence note.
  `other_notes`: any other statement on a witness's presence (beginning, end, lacuna, fragment), verbatim.
- One unit per lemma. Witness sigla as printed (`F14`, `L1`, `U1`, `SO`, `Gb`, `Nn`); printed editions as `BUON`,
  `DRAG`, `LAZZ`, `MOMB`.
- Group sigla are written as ASCII names (`kappa`, `lambda`, `mu`, `nu`, `xi`, `omicron`, `pi`, `rho`, `sigma`,
  `tau`, `upsilon`, `phi`), with exceptions as `"nu(-rho,-V3)"` for "ν (tranne ϱ, V3)". The builder expands each
  group into its members.
- Corrections are encoded by hand: `Li^ac` (first hand) under the reading before correction, `Li^pc` under the
  reading after it; additions above the line and in the margin likewise.
- `type` is one of `spelling`, `word`, `order`, `omission`, `addition`, `other`; `uncertain` marks any siglum,
  digit or group letter not certain on the page.
