"""Vision-LLM page reader (A8 candidate): full-page verbatim transcription of an MGH page image, line by line.
usage: bedier_llm_read.py MODEL VOL LEAF OUTDIR   (A8 source; prompt and settings frozen)"""
import os, sys, json, base64, urllib.request, time
model, vol, leaf, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
key = os.environ['OPENROUTER_API_KEY']
imgd = os.environ.get('DMGH_IMG', os.path.join(out, 'img')); os.makedirs(imgd, exist_ok=True)
ip = os.path.join(imgd, f'{vol}_{leaf:05d}.jpg')
if not os.path.exists(ip):
    urllib.request.urlretrieve(f'https://www.dmgh.de/{vol}/img/300/{vol}_{leaf:05d}.jpg', ip)
PROMPT = """This is a page of a 19th-century Monumenta Germaniae Historica edition (Latin). Transcribe the WHOLE page exactly as printed, line by line, top to bottom: running head with page number, edited text with its marginal line numbers, then the critical apparatus at the foot (small type).
Rules:
- One printed line per output line; keep line numbers, brackets ']', punctuation, parentheses, italics as plain text.
- The apparatus is critical: copy every manuscript siglum exactly (capital letters, with superscript hands written as plain digits after the letter, e.g. P2, C1; lower-case superscripts as plain lower-case, e.g. Mb). Keep runs of sigla exactly as printed (e.g. LMTCFP). Never add, normalise or guess a siglum; if illegible write ?.
- If the apparatus is printed in two columns, transcribe the left column fully, then the right column.
- Output only the transcription, no comments, no markdown."""
b64 = base64.b64encode(open(ip, 'rb').read()).decode()
body = {'model': model, 'temperature': 0, 'messages': [{'role': 'user', 'content': [
    {'type': 'text', 'text': PROMPT}, {'type': 'image_url', 'image_url': {'url': 'data:image/jpeg;base64,' + b64}}]}]}
for attempt in range(4):
    try:
        req = urllib.request.Request('https://openrouter.ai/api/v1/chat/completions', data=json.dumps(body).encode(),
                                     headers={'Authorization': 'Bearer ' + key, 'Content-Type': 'application/json'})
        r = json.load(urllib.request.urlopen(req, timeout=600)); break
    except Exception as e:
        time.sleep(5 * (attempt + 1)); err = e
else: raise err
t = r['choices'][0]['message']['content'].strip().strip('`'); u = r.get('usage', {})
os.makedirs(os.path.join(out, vol), exist_ok=True)
open(os.path.join(out, vol, f'{vol}_{leaf:05d}.llm.txt'), 'w').write(t + '\n')
with open(os.path.join(out, 'USAGE.tsv'), 'a') as fh: fh.write(f"{vol}\t{leaf}\t{r.get('model')}\t{u.get('prompt_tokens')}\t{u.get('completion_tokens')}\t{u.get('cost')}\n")
print(vol, leaf, r.get('model'), u.get('prompt_tokens'), u.get('completion_tokens'), u.get('cost'))
