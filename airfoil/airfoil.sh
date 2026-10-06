#!/usr/bin/env bash
# Usage: ./run_polars.sh            (runs every .dat in the current dir)
#        RES="200000 500000" ./run_polars.sh
set -u

OUT_DIR="${OUT_DIR:-polars}"
RES="${RES:-500000 1000000}"
A_MAX="${A_MAX:-15}"
A_MIN="${A_MIN:--6}"
A_STEP="${A_STEP:-1}"
NCRIT="${NCRIT:-9}"
ITER="${ITER:-200}"
TIMEOUT="${TIMEOUT:-7}"     # hard wall-clock limit per run (s)

mkdir -p "$OUT_DIR"
FAILED="$OUT_DIR/failed.txt"; : > "$FAILED"

# Never leave a stray xfoil behind on exit or Ctrl-C.
# (Remove the pkill if you run several copies of this script at once.)
tmp=""
cleanup() { pkill -9 -x xfoil 2>/dev/null; [ -n "$tmp" ] && rm -rf "$tmp"; }
trap 'cleanup; exit 130' INT TERM
trap cleanup EXIT

pkill -9 -x xfoil 2>/dev/null    # clear anything left over from earlier runs

for f in ./*.dat; do
  name=$(basename "$f" .dat)
  echo "$name"
  for Re in $RES; do
    out="$OUT_DIR/${name}_Re${Re}.txt"
    [ -s "$out" ] && continue

    tmp=$(mktemp -d)
    cp "$f" "$tmp/$name.dat"

    cat > "$tmp/cmds.in" <<EOF

load ${name}.dat
pane
oper
visc $Re
mach 0
vpar
n $NCRIT

iter $ITER
pacc
polar.txt

aseq 0 $A_MAX $A_STEP
init
aseq -$A_STEP $A_MIN -$A_STEP
pacc

quit
quit
quit
quit
EOF

    # SIGTERM at $TIMEOUT, then SIGKILL 3s later if it's still alive.
    ( cd "$tmp" && timeout -k 3 "$TIMEOUT" xfoil < cmds.in > xfoil.log 2>&1 )
    status=$?
    pkill -9 -x xfoil 2>/dev/null

    if [ -s "$tmp/polar.txt" ] && [ "$(wc -l < "$tmp/polar.txt")" -gt 12 ]; then
      { head -n 12 "$tmp/polar.txt"; tail -n +13 "$tmp/polar.txt" | sort -g -k1,1; } > "$out"
      if [ "$status" -eq 124 ] || [ "$status" -eq 137 ]; then
        echo "OK*   $name Re=$Re (killed after timeout, partial polar kept)"
        echo "$name Re=$Re PARTIAL exit=$status" >> "$FAILED"
      else
        echo "OK    $name Re=$Re"
      fi
    else
      echo "$name Re=$Re exit=$status" >> "$FAILED"
      cp "$tmp/xfoil.log" "$OUT_DIR/${name}_Re${Re}.log"
      echo "FAIL  $name Re=$Re (exit $status; 124/137 = timeout)"
    fi
    rm -rf "$tmp"; tmp=""
  done
done

echo "Done. Failures listed in $FAILED"