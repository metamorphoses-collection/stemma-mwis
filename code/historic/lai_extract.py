"""Lai de l'Ombre, Bedier 1913 (SATF): extract shared departures from the vision transcriptions of the variant notes
(data/historic/lai_ombre/vision/nNNN.txt). Entry format: '<verse> SIGLA reading, SIGLA reading — <verse> ...'. A run of
sigla A-G opening a reading (at the start of an entry or after a comma) is a group of witnesses sharing a departure
from the printed text (manuscript A, corrected in places). Parenthesised sub-readings are dropped (as in the MGH parser);
'omet' readings count as departures. Output: units.json in the format read by bedier_top.py."""
import re, json, sys, glob, os, collections
SIG = r'[A-G]{1,7}'


def entries(text):
    text = re.sub(r'\([^()]*\)', ' ', text)
    for ent in re.split(r'\s+—\s+|\s—|—\s', text):
        m = re.match(r'^\s*(\d{1,3})(?:\s*-\s*(\d{1,3}))?\s+(.*)$', ent.strip())
        if not m: continue
        verse, body = int(m.group(1)), m.group(3)
        groups = []
        for g in re.finditer(r'(?:^|,\s+|;\s+)(' + SIG + r')(?=\s|$)', body):
            s = g.group(1)
            if len(set(s)) == len(s): groups.append(''.join(sorted(s)))
        if groups: yield verse, groups, ent.strip()[:200]


def main(vdir, out):
    units = []
    for f in sorted(glob.glob(os.path.join(vdir, 'n*.txt'))):
        t = open(f).read()
        if 'NOTES:' not in t or 'NOTES: none' in t: continue
        notes = t.split('NOTES:', 1)[1].replace('\n', ' ')
        for verse, groups, raw in entries(notes):
            # same witness twice in one entry (different parts of the verse): keep each group, as separate readings
            units.append(dict(page=os.path.basename(f)[:-4], line=verse, entry_line=verse, leaf=None, lemma=None, lemma_group=0,
                              groups=[''] + groups, readings=[raw]))
    json.dump(units, open(out, 'w'), ensure_ascii=False, indent=0)
    c = collections.Counter(g for u in units for g in u['groups'][1:] if len(g) >= 2)
    print(len(units), 'entries;', sum(c.values()), 'shared departures in', len(c), 'groups')
    print(c.most_common(30))
    print("Bedier's checks: DEF", c['DEF'], '+ ABCG', c['ABCG'], '=', c['DEF'] + c['ABCG'], '(he says ~30) | CDEFG', c['CDEFG'], '+ AB', c['AB'], '=', c['CDEFG'] + c['AB'], '(he says ~80)')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
