#!/bin/bash
# Simple Tor-based GP downloader with circuit rotation
# Usage: ./fetch-tor.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../tests/data/socrates_gp"
CSV="$SCRIPT_DIR/../tests/data/socrates_current.csv"

mkdir -p "$DATA_DIR"

# Extract pairs not yet downloaded
echo "Finding missing pairs..."
tail -n +2 "$CSV" | awk -F',' '{gsub(/"/, ""); print $1","$4}' | sort -u > /tmp/all_pairs.txt
TOTAL=$(wc -l < /tmp/all_pairs.txt | tr -d ' ')
echo "Total unique pairs: $TOTAL"

# Filter out already downloaded
> /tmp/need_pairs.txt
while IFS=',' read -r id1 id2; do
    [ -f "$DATA_DIR/gp_${id1},${id2}.json" ] || echo "${id1},${id2}" >> /tmp/need_pairs.txt
done < /tmp/all_pairs.txt
NEED=$(wc -l < /tmp/need_pairs.txt | tr -d ' ')
echo "Already have: $((TOTAL - NEED))"
echo "Need to download: $NEED"
echo ""

rotate() {
    (echo -e 'AUTHENTICATE ""\r\nSIGNAL NEWNYM\r\nQUIT\r') | nc -w 2 127.0.0.1 9051 >/dev/null 2>&1 || true
    sleep 2
}

dl=0; err=0; rot=0; cerr=0

while IFS=',' read -r id1 id2; do
    out="$DATA_DIR/gp_${id1},${id2}.json"
    [ -f "$out" ] && continue

    code=$(curl --socks5 127.0.0.1:9050 -s -w "%{http_code}" -o "$out.tmp" \
        --max-time 12 -H "User-Agent: Mozilla/5.0" \
        "https://celestrak.org/SOCRATES/data.php?CATNR=${id1},${id2}&FORMAT=json" 2>/dev/null || echo "000")

    if [ "$code" = "200" ] && python3 -c "import json; d=json.load(open('$out.tmp')); assert len(d)>=2" 2>/dev/null; then
        mv "$out.tmp" "$out"
        dl=$((dl+1)); cerr=0
    elif [ "$code" = "429" ] || [ "$code" = "503" ]; then
        rm -f "$out.tmp"; rot=$((rot+1))
        echo "  ⟳ 429/503 at $dl+$err, rotating ($rot)"
        rotate; sleep 2; continue
    else
        rm -f "$out.tmp"; err=$((err+1)); cerr=$((cerr+1))
        if [ $cerr -ge 8 ]; then
            rot=$((rot+1)); echo "  ⟳ streak=$cerr, rotating ($rot)"
            rotate; cerr=0
        fi
    fi

    tot=$((dl+err))
    [ $((tot % 500)) -eq 0 ] && [ $tot -gt 0 ] && echo "  $dl downloaded, $err errors, $rot rotations (of $NEED)"
    sleep 0.06
done < /tmp/need_pairs.txt

echo ""
echo "Done: $dl downloaded, $err errors, $rot rotations"
echo "Total GP files: $(ls "$DATA_DIR"/*.json | wc -l)"
