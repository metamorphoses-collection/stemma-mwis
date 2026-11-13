# Verification of the recorded editors' top divisions (before freezing)
Each agent reading is checked by the lead session against the quoted Latin at the given OCR line. Status: accepted /
amended (with reason) / query. Readings are fixed when this file and tops/*.json are committed; that commit hash is the
one cited for the blind step.

| Volume | Text | Recorded top | Check | Status |
|---|---|---|---|---|
| AA 1 | Salvianus, De gub. Dei / Ad eccl. | none stated | agent: Halm states no grouping | accepted (excluded: no stemma) |
| AA 1 | Eugippius, Vita Severini (Sauppe) | L / VMD | p. XV l.15119 "L et VM ex duobus apographis unius archetypi" ✓; D with VM rests on p. X l.14867 "magis cum codicibus VM quam cum L convenire" (reader inference, D a 4-leaf fragment, likely below 30 citations); aside p. XIII l.15003 on an LM error noted | accepted, D marked inferred |
| AA 3 | Victor Vitensis, Historia (Halm) | BLV / R Bern. p Lor. | p. V l.62 "genera sunt duo" ✓; p. VI l.112 BLV ✓; p. VII l.147 (Rp) ✓; A (Laudunensis) unplaced, p. IX l.252 "modo in BV, modo in R" | accepted; A outside the top branches |
| AA 3 | Corippus, In laudem Iustini (Partsch) | M / (O S L excerpts) | diagram p. LIX garbled in OCR; covers III 271–407 only; O S L indirect | query: check diagram on image; likely excluded (too few direct witnesses) |
| AA 7 | Ennodius (Vogel) | B / VLTCOP | p. XXXII "ab altera parte stat unus codex Bruxellensis, ab altera ceteri omnes" (agent quote; repeated p. XLV) | accepted; diagram p. XLVIII to view on image |
| AA 8 | Sidonius, letters (Luetjohann/Leo) — PILOT | L / MTCFP | stemma p. XLI (α → L, β); β completed from σ in gap passages (flag contamination_at_top) | accepted (pilot, reported separately) |
| AA 8 | Sidonius, carmina (Leo) | M / TCFP | σ (M, carm. I–VIII only) vs x; T's source for I–VIII left open by Leo | accepted, flagged; top visible only in carm. I–VIII |
| AA 8 | Faustus et al., letters (Krusch) | none stated | no filiation | accepted (excluded: no stemma) |
| AA 12 | Cassiodorus, Variae — PILOT | MONXZEF / LRKPD (I–IV,39) | classes by books from one archetype; two classes side by side only in I–IV,39 | accepted (pilot) |
| AA 12 | Acta synhodorum 499/501/502 | ABFT / EC | p. 395 l.49020 "exemplaria duo … A … E" ✓; p. 397 l.49173 third order "contaminatorum" (H G D M) outside the two branches; K V excerpts from a different codex (indirect) | accepted; HGDM unplaced (contaminated) |
| AA 13 | Gildas, De excidio | C(D) / AX | p. 13 l.1114 all "redeunt fere ad unum codicem … Cottonianum" except interpolated; p. 17 l.1393 "archetypum interpolatum" ✓; D descriptus of C | accepted |
| AA 13 | Historia Brittonum | HKPQDGCL / ZMN | p. 117 l.10693 "primae tertiaeque classis archetypo communi" ✓ → I+III vs II; p. 135 l.12058 heads of II/III "nequaquam descendere" from I's archetype; MN contaminated from the Gildan order (p. 136) | accepted as two-way (editor's words); flagged contamination |
| AA 13 | Beda, Chronica maiora / minora | none | p. 243–244 "familiae codicum ordinari nequiverunt" | accepted (excluded: no stemma) |
| AA 13 | Laterculi (imperatorum; regum Visigothorum) | none | only local groups | accepted (excluded: no stemma) |
| AA 4,1 | Fortunatus, Vita Martini (Leo) | **amended** PTF / NBLMS | printed stemma p. XXIII read on image (leaf 22): Fort. → α (P; x→T,F), β (N; γ interp.→BL, MS). Agent had P / NBLMS (T, F unplaced) | amended to the diagram |
| AA 4,2 | Vita Radegundis (Krusch) | PMG / LEW / TH+XYZ (k=3) | image p. XVIII (leaf 17): "Interpolatos in tres classes distribui has: LEW — TH — XYZ"; "Libri TH cum LE vinculo nullo coniuncti sunt, sed XY procul dubio e TH fluxerunt" | accepted, three-branch top |
| AA 4,2 | Vita Medardi (Krusch) | MNQB / PR | p. XXVIII "Codicum duo genera distinguuntur, quorum alteri libros MNQB, alteri PR tribui" (seen on leaf 27 OCR) | accepted |
| AA 4,2 | Vita Hilarii | GHB / VA | H from G, B from H; A "praecedenti [V] persimilis"; two branches implied, not stated | query: implied only |
| AA 4,2 | Vita Germani | PQ / GHI | "Q cum P semper convenit … HI et G2 consentiunt"; KL contaminated from P and G2 | query: implied only |
| AA 4,2 | Expositio fidei | A / MF Oxon. | MF share errors; Oxon. "eiusdem generis ac MF" | query: implied only |
| AA 4,2 | Vitae Albini, Paterni, Marcelli, Amantii, Remedii, Leobini, Maurilii; Passio Dionysii | none | no filiation stated | accepted (excluded: no stemma) |
| AA 4,1 | Fortunatus, carmina — PILOT | Σ / rest | all except Σ from one mutilated archetype (p. XVII); Σ a selection | accepted (pilot) |
| AA 2 | Eutropius — PILOT | GFD(P) / LO(Paean.) | A+C archetype vs B (p. XIII l.1332) | accepted (pilot) |
| AA 2 | Paulus, Historia Romana | HVBM / NZ | two kinds distinguished by the lacuna in book XV (p. XXIX); M corrected from 2nd family | accepted, flag contamination |
| AA 14 | Eugenius, carmina (Vollmer) | **amended** Λ / M / L / Brux. / Harl.+Lips. / F P C R T (k=6) | p. XLV "nihil conicere licet de superiore singulorum affinitate"; printed stemma p. XLVI read on image (leaf 47): Eug. → Λ, M, L, Brux., Π(Harl., Lips.), Φ(F, Ξ(P, Y(Z(C,R), T))) | amended to the diagram |
| AA 10 | Claudian, carmina maiora (Birt) | GVP(M) / rest | stemma p. CIII on image (leaf 106): Ω → x (G, V, P) | y; P also from w | accepted, P flagged contaminated; M to recheck |
| AA 10 | Claudian, De raptu Proserpinae | FVS / CABDW | p. CXLVIII two recensions; p. CXLIX "media" class "cum altera nexum artissimum" | accepted (two-way) |
| AA 10 | Claudian, Paneg. Probino et Olybrio | none | only a ranking BRW better than LFPT | accepted (excluded: no stemma) |
| AA 14 | Dracontius (Laudes, Satisfactio, Orestis) | as recorded | B only full direct witness / 2–3 direct witnesses | accepted; will fail the witness rule |
| AA 14 | Eugenius, Dracontiana | M / FPLQZ | footnote anchor garbled | query: check pp. XIX–XXII on image |
| AA 15 | Aldhelm (Ehwald), five texts | as recorded | two classes stated for De metris, Aenigmata (recensions), De virginitate prosa and carmen; Carmina eccl.: S "suo stat talo" | accepted; Carmina eccl. implied only |
| AA 11 | Chronica minora II: Marcellinus; Cassiodorus chron.; Victor; Iohannes Biclar.; Isidore Hist. Goth.; Isidore Chron. maiora | as recorded | Isidore: two recensions/editions; Marcellinus, Cassiodorus: two independent witnesses only | accepted; most fail the witness rule |
| AA 9 | Prosper, Epitoma chronicon | AORV / MYL | p. 376 two primary recensions, "neutra ab altera pendeat"; CDHBFPXZ contaminated from both (flag) | accepted; contaminated witnesses outside the top branches |
| AA 9 | Consularia Constantinopolitana | C / Sm Sv Sp | S group = three copies of one lost codex | accepted; fails witness rule (2 independent) |
| AA 9 | Chronica Gallica 452 | L / BM | p. 620 two "ordines" | accepted; fails witness rule (3) |
| AA 9 | Liber generationis I | BFO / GC | pp. 78–81 | accepted |
| AA 9 | Chronographus 354 | star, "potest defendi" | tentative derivation from the lost Luxemburg codex | **undecided** (reported separately) |
| AA 9 | Computatio 452 | none | S "fortasse" archetype of all | accepted (excluded: no top division) |
| AA 9 | Victorius, Cursus paschalis | — | a ranking with shared interpolations, not an archetype split; C wavers | **undecided** |
| AA 9 | Polemius, Laterculus | — | branch placement is the agent's reading of p. 524 | **undecided** unless confirmed on image |
| AA 9 | Notitia Galliarum; Liber genealogus | none | p. 560 no true families because of contamination; genealogus "non recensui, sed repraesentavi" | accepted (excluded: no stemma) |
| AA 5 | Jordanes, Getica — PILOT | HPVLAOB(S) / XYZ | p. LXXII; l.4005 "ab una parte stant HPVLASOB contra XYZ" | accepted (pilot) |
| AA 5 | Jordanes, Romana | HPVL / S (indirect) | class II only in excerpts | accepted; one direct branch, nothing to classify |
| AA 5 | Ausonius (Schenkl) | — | "duo sunt codicum genera", but each a set of excerpts; several archetypes | **undecided** |
| AA 6,1 | Symmachus, letters (Seeck) | P(Π) / VMF ΓΦ | full stemma p. XXXVII on image (leaf 36): Archetypus → P (→ Π) | rest (P2m → V, M; F → F1–3; Γ; Φ); also p. XXXII (leaf 31) | accepted |
| AA 6,1 | Symmachus, relationes | TM / Γ (Gelenius ed., indirect) | p. XX | accepted; one direct branch |
| AA 6,2 | Avitus, letters and homilies (Peiper) | papyrus / L P S | stemma p. XXVIII on image (leaf 27): three branches (Gregory's exemplar, lost; Lugdunensis → L, PS; papyrus) | accepted |
| AA 6,2 | Avitus, carmina | — | p. LII on image (leaf 51): three ordines within the Gallicani only; relation to archetype "incertissimum"; no α/β top stated | **undecided** |
