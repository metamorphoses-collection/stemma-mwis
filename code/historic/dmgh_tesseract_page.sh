#!/bin/bash
# one.sh VOL LEAF OUTDIR : download one dMGH leaf at 300 dpi, OCR with Tesseract (Latin), keep text only
v=$1; l=$(printf "%05d" $2); o=$3
[ -s "$o/$v/$l.txt" ] && exit 0
mkdir -p "$o/$v"
img="$o/$v/$l.jpg"
for t in 1 2 3; do curl -s --max-time 90 -o "$img" "https://www.dmgh.de/$v/img/300/${v}_$l.jpg" && [ -s "$img" ] && break; sleep 3; done
tesseract "$img" "$o/$v/$l" -l lat --psm 6 >/dev/null 2>&1
rm -f "$img"
