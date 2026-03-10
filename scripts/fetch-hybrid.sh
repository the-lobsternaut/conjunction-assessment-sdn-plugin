#!/bin/bash
# Hybrid GP downloader: direct first, Tor on rate limit
# Usage: ./fetch-hybrid.sh [parallel]

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../tests/data/socrates_gp"
CSV="$SCRIPT_DIR/../tests/data/socrates_current.csv"
PARALLEL="${1:-8}"
USE_TOR=0
DL=0
ERR=0
ROT=0

mkdir -p "$DATA_DIR"

# Build needed list
python3 - "$DATA_DIR" "$CSV" << 'PYEOF'
import os, csv, sys
data_dir, csv_path = sys.argv[1], sys.argv[2]
existing = set(f.replace('gp_', '').replace('.json', '') for f in os.listdir(data_dir) if f.endswith('.json'))
seen = set()
needed = []
with open(csv_path) as f:
    for row in csv.DictReader(f):
        key = f"{row['NORAD_CAT_ID_1']},{row['NORAD_CAT_ID_2']}"
        if key not in seen and key not in existing:
            seen.add(key)
            needed.append(key)
with open('/tmp/gp_needed.txt', 'w') as f:
    for p in needed: f.write(p + '\n')
print(f"Have: {len(existing)}, Need: {len(needed)}")
PYEOF

NEED=$(wc -l < /tmp/gp_needed.txt | tr -d ' ')
echo "Downloading $NEED GP pairs (${PARALLEL} parallel, Tor fallback)..."

fetch_pair() {
    local pair="$1"
    local id1="${pair%%,*}"
    local id2="${pair#*,}"
    local out="$DATA_DIR/gp_${id1},${id2}.json"
    [ -f "$out" ] && return 0

    local proxy_flag=""
    [ -f /tmp/use_tor ] && proxy_flag="--socks5 127.0.0.1:9050"

    local code
    code=$(curl $proxy_flag -s -w "%{http_code}" -o "$out.tmp" --max-time 10 \
        -H "User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36" \
        "https://celestrak.org/SOCRATES/data.php?CATNR=${id1},${id2}&FORMAT=json" 2>/dev/null || echo "000")

    if [ "$code" = "200" ] && python3 -c "import json; d=json.load(open('$out.tmp')); assert len(d)>=2" 2>/dev/null; then
        mv "$out.tmp" "$out"
        return 0
    fi
    rm -f "$out.tmp"

    if [ "$code" = "429" ] || [ "$code" = "503" ]; then
        # Signal to switch to Tor
        touch /tmp/use_tor
        (echo -e 'AUTHENTICATE ""\r\nSIGNAL NEWNYM\r\nQUIT\r') | nc -w 2 127.0.0.1 9051 >/dev/null 2>&1 || true
        sleep 2
        
        # Retry via Tor
        code=$(curl --socks5 127.0.0.1:9050 -s -w "%{http_code}" -o "$out.tmp" --max-time 15 \
            -H "User-Agent: Mozilla/5.0" \
            "https://celestrak.org/SOCRATES/data.php?CATNR=${id1},${id2}&FORMAT=json" 2>/dev/null || echo "000")
        
        if [ "$code" = "200" ] && python3 -c "import json; d=json.load(open('$out.tmp')); assert len(d)>=2" 2>/dev/null; then
            mv "$out.tmp" "$out"
            return 0
        fi
        rm -f "$out.tmp"
    fi
    return 1
}
export -f fetch_pair
export DATA_DIR

rm -f /tmp/use_tor

# Run parallel downloads
cat /tmp/gp_needed.txt | xargs -P "$PARALLEL" -I {} bash -c 'fetch_pair "$@"' _ {} &
BG=$!

# Monitor
while kill -0 $BG 2>/dev/null; do
    sleep 15
    count=$(find "$DATA_DIR" -name "*.json" -maxdepth 1 | wc -l | tr -d ' ')
    mode="direct"
    [ -f /tmp/use_tor ] && mode="TOR"
    echo "  [$mode] $count GP files"
done

wait $BG 2>/dev/null
FINAL=$(find "$DATA_DIR" -name "*.json" -maxdepth 1 | wc -l | tr -d ' ')
echo ""
echo "Done! Total: $FINAL GP files"
rm -f /tmp/use_tor
