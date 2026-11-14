"""Vision transcription of one page of Bedier's Lai de l'Ombre (SATF 1913, Internet Archive images): the variant notes
at the foot of the page, verbatim. usage: lai_read.py N OUTDIR  (N = Internet Archive page index)"""
import os, sys, json, base64, urllib.request, time
n, out = int(sys.argv[1]), sys.argv[2]
os.makedirs(out, exist_ok=True); img = os.path.join(out, f'n{n:03d}.jpg'); txt = os.path.join(out, f'n{n:03d}.txt')
if os.path.exists(txt): sys.exit()
if not os.path.exists(img): urllib.request.urlretrieve(f'https://archive.org/download/lelaidelombre00jeanuoft/page/n{n}.jpg', img)
P = ("This is a page of Joseph Bedier's 1913 edition of the Old French Lai de l'Ombre. Transcribe exactly: first line the "
     "running head with the page number, then the line 'VERSES:' followed by the verse numbers printed in the margin "
     "(first and last on the page, if any), then the line 'NOTES:' followed by the variant notes at the foot of the page, "
     "verbatim, keeping every verse number, every dash (—) and every manuscript siglum. The sigla are the capital letters "
     "A B C D E F G, printed in italics, and they PRECEDE the reading they introduce (e.g. '447 CG K'il saiche (G Qui sace)'). "
     "Keep runs of sigla exactly as printed (CDEFG). Never add or normalise a siglum. If there are no notes, write 'NOTES: none'. "
     "Output only the transcription.")
body = {'model': 'google/gemini-4-argon', 'temperature': 0, 'messages': [{'role': 'user', 'content': [
    {'type': 'text', 'text': P}, {'type': 'image_url', 'image_url': {'url': 'data:image/jpeg;base64,' + base64.b64encode(open(img, 'rb').read()).decode()}}]}]}
for a in range(4):
    try:
        req = urllib.request.Request('https://openrouter.ai/api/v1/chat/completions', data=json.dumps(body).encode(),
                                     headers={'Authorization': 'Bearer ' + os.environ['OPENROUTER_API_KEY'], 'Content-Type': 'application/json'})
        r = json.load(urllib.request.urlopen(req, timeout=600)); break
    except Exception as e: err = e; time.sleep(5 * (a + 1))
else: raise err
open(txt, 'w').write((r['choices'][0]['message']['content'] or '').strip() + '\n')
