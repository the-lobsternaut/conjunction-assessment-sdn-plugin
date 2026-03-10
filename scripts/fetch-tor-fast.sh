#!/bin/bash
# Fast Tor-based GP downloader — parallel with xargs, circuit rotation on errors
# Usage: ./fetch-tor-fast.sh [parallel_jobs]

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../tests/data/socrates_gp"
CSV="$SCRIPT_DIR/../tests/data/socrates_current.csv"
JOBS="${1:-4}"

mkdir -p "$DATA_DIR"

# Build list of needed pairs quickly using python
echo "Building download list..."
python3 << 'PYEOF'
import os, csv

data_dir = os.environ.get('DATA_DIR', 'tests/data/socrates_gp')
csv_path = os.environ.get('CSV', 'tests/data/socrates_current.csv')

existing = set(f.replace('gp_', '').replace('.json', '') for f in os.listdir(data_dir) if f.endswith('.json'))
print(f"Existing: {len(existing)}", flush=True)

seen = set()
needed = []
with open(csv_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = f"{row['NORAD_CAT_ID_1']},{row['NORAD_CAT_ID_2']}"
        if key not in seen and key not in existing:
            seen.add(key)
            needed.append(key)

print(f"Need: {len(needed)}", flush=True)
with open('/tmp/gp_needed.txt', 'w') as f:
    for pair in needed:
        f.write(pair + '\n')
PYEOF

NEED=$(wc -l < /tmp/gp_needed.txt | tr -d ' ')
echo "Downloading $NEED pairs via Tor (${JOBS} parallel)..."

# Download function
download_one() {
    pair="$1"
    id1="${pair%%,*}"
    id2="${pair#*,}"
    out="$DATA_DIR/gp_${id1},${id2}.json"
    [ -f "$out" ] && return 0

    for attempt in 1 2 3; do
        code=$(curl --socks5 127.0.0.1:9050 -s -w "%{http_code}" -o "$out.tmp" \
            --max-time 15 -H "User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)" \
            "https://celestrak.org/SOCRATES/data.php?CATNR=${id1},${id2}&FORMAT=json" 2>/dev/null || echo "000")

        if [ "$code" = "200" ]; then
            if python3 -c "import json; d=json.load(open('$out.tmp')); assert len(d)>=2" 2>/dev/null; then
                mv "$out.tmp" "$out"
                return 0
            fi
        fi
        rm -f "$out.tmp"

        if [ "$code" = "429" ] || [ "$code" = "503" ]; then
            # Rate limited — rotate and retry
            (echo -e 'AUTHENTICATE ""\r\nSIGNAL NEWNYM\r\nQUIT\r') | nc -w 2 127.0.0.1 9051 >/dev/null 2>&1 || true
            sleep $((attempt * 2))
        else
            sleep 0.5
        fi
    done
    return 1
}
export -f download_one
export DATA_DIR

# Run with xargs for parallelism, monitoring progress
cat /tmp/gp_needed.txt | xargs -P "$JOBS" -I {} bash -c 'download_one "$@"' _ {} &
XARGS_PID=$!

# Monitor progress
while kill -0 $XARGS_PID 2>/dev/null; do
    sleep 30
    count=$(ls "$DATA_DIR"/*.json 2>/dev/null | wc -l | tr -d ' ')
    echo "  Progress: $count GP files"
done

wait $XARGS_PID 2>/dev/null
FINAL=$(ls "$DATA_DIR"/*.json 2>/dev/null | wc -l | tr -d ' ')
echo ""
echo "Done! Total GP files: $FINAL"
