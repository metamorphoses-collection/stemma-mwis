# Bédier's own case: the Lai de l'Ombre (exploratory, 14 Nov 2026)
Source: Bédier's 1913 edition (SATF), Internet Archive `lelaidelombre00jeanuoft`, pages n55–n115 transcribed by the
vision model (code/historic/lai_read.py; data/historic/lai_ombre/vision/). Parser: code/historic/lai_extract.py
(sigla precede readings; departures from the printed text, which is manuscript A corrected). 814 notes, verses 1–957.
Check against Bédier's own tallies: DEF vs ABCG 24 (he: ~30); identical CDEFG vs AB 38 (he: ~80 oppositions, counted
whether or not C, G, D, E, F agree among themselves, so a lower identical count is expected).
Rooted run (bedier_top.py, B = 200): stable CG 176 (1.00), DE 43 (0.995), CDEFG 21 (1.00), DEF 21 (0.995), AB 17 (1.00).
- Bédier 1890 (Schéma 1, ABCG | DEF): contradicted by CDEFG (rooted by the A text).
- Paris 1890 (Schéma 2, ABCG | DF | E): contradicted by CDEFG and DE.
Unrooted split weights: CG 176, DE 43, AB|CDEFG 38, ABCG|DEF 24; EF 23 and DF 13 cross DE. Tree ((A,B),(C,G),((D,E),F)).
Bédier 1890, Quentin (AB | CGDEF) and a three-branch top AB | CG | DEF are this tree rooted at three places; the two
candidate root edges weigh 21 each against the printed text. The readings fix the groups; the root needs a judgement on
errors. Bédier (1928, pp. 339–351) reached the same conclusion: schémas 1, 3, 8–11 are all "non moins plausibles".
Sources: Bédier 1928 (Romania 54), Quentin 1926, Paris 1890; Roelli (ed.) 2020, sect. 2.3 and ex. 6 of 6.2.

## Proof page (blind)
Drawn with random.Random("LAI-proof-20261114") among the 47 pages with ≥ 8 notes, excluding n080 (seen before):
n103 = printed p. 46. Gold transcribed blind from the image by an AI agent (proof_gold_n103.json).
Shared groups (≥ 2 sigla): 4 printed, 4 recovered, 0 extra. All readings incl. single witnesses, by verse: 42 of 43,
0 extra; the miss is v. 914 "vers omis en G" (siglum after the text), a single-witness omission. proof_score_n103.json.
