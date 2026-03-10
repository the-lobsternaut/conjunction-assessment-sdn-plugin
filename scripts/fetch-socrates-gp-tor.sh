#!/bin/bash
# Batch download SOCRATES GP data via Tor with circuit rotation
# Usage: ./fetch-socrates-gp-tor.sh [max_pairs] [rate_ms]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../tests/data/socrates_gp"
CSV_FILE="$SCRIPT_DIR/../tests/data/socrates_current.csv"
MAX_PAIRS="${1:-50000}"
RATE_SEC="${2:-0.08}"  # 80ms between requests = ~12/sec

mkdir -p "$DATA_DIR"

rotate_circuit() {
    # Send NEWNYM signal to Tor control port
    (echo -e 'AUTHENTICATE ""\r\nSIGNAL NEWNYM\r\nQUIT\r') | nc -w 3 127.0.0.1 9051 >/dev/null 2>&1 || true
    sleep 2  # Wait for new circuit
}

# Extract unique pairs from CSV
echo "Parsing SOCRATES CSV..."
PAIRS=$(tail -n +2 "$CSV_FILE" | head -n "$MAX_PAIRS" | awk -F',' '{
    # Handle quoted fields
    gsub(/"/, "", $0)
    id1=$1; id2=$4
    key=id1","id2
    if (!(key in seen)) {
        seen[key]=1
        print id1","id2
    }
}')

TOTAL=$(echo "$PAIRS" | wc -l | tr -d ' ')
echo "Total unique pairs: $TOTAL"
echo "Output: $DATA_DIR"
echo "Rate: ${RATE_SEC}s between requests"
echo ""

downloaded=0
cached=0
errors=0
rotations=0
consecutive_errors=0
count=0

while IFS=',' read -r id1 id2; do
    count=$((count + 1))
    outfile="$DATA_DIR/gp_${id1},${id2}.json"
    
    # Skip if already downloaded
    if [ -f "$outfile" ]; then
        cached=$((cached + 1))
        continue
    fi
    
    url="https://celestrak.org/SOCRATES/data.php?CATNR=${id1},${id2}&FORMAT=json"
    
    # Download via Tor
    http_code=$(curl --socks5 127.0.0.1:9050 -s -w "%{http_code}" \
        -o "$outfile.tmp" --max-time 15 \
        -H "User-Agent: Mozilla/5.0 (compatible)" \
        "$url" 2>/dev/null || echo "000")
    
    if [ "$http_code" = "200" ]; then
        # Verify it's valid JSON with 2+ objects
        if python3 -c "import json; d=json.load(open('$outfile.tmp')); assert len(d) >= 2" 2>/dev/null; then
            mv "$outfile.tmp" "$outfile"
            downloaded=$((downloaded + 1))
            consecutive_errors=0
        else
            rm -f "$outfile.tmp"
            errors=$((errors + 1))
            consecutive_errors=$((consecutive_errors + 1))
        fi
    elif [ "$http_code" = "429" ] || [ "$http_code" = "503" ]; then
        rm -f "$outfile.tmp"
        rotations=$((rotations + 1))
        echo "  ⟳ Rate limited at $count ($id1,$id2), rotating circuit... ($rotations)"
        rotate_circuit
        sleep 3
        # Retry — decrement count so we try again
        count=$((count - 1))
        continue
    else
        rm -f "$outfile.tmp"
        errors=$((errors + 1))
        consecutive_errors=$((consecutive_errors + 1))
    fi
    
    # Rotate on error streaks
    if [ "$consecutive_errors" -ge 10 ]; then
        rotations=$((rotations + 1))
        echo "  ⟳ $consecutive_errors consecutive errors, rotating... ($rotations)"
        rotate_circuit
        consecutive_errors=0
    fi
    
    # Progress
    processed=$((downloaded + errors))
    if [ $((processed % 500)) -eq 0 ] && [ "$processed" -gt 0 ]; then
        echo "  [$count/$TOTAL] $downloaded new, $cached cached, $errors err, $rotations rotations"
    fi
    
    sleep "$RATE_SEC"
done <<< "$PAIRS"

echo ""
echo "Done: $downloaded downloaded, $cached cached, $errors errors, $rotations rotations"
echo "Total GP files: $((downloaded + cached))"
