#!/bin/bash

ARGS="$@"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

LOGDIR="DAP_Logs"
LOGFILE="$LOGDIR/DAP_$1_$2_${TIMESTAMP}.log"

mkdir -p "$LOGDIR"

python fit_one_cube_muse.py "$ARGS" > "$LOGFILE" 2>&1 &

echo "Logging output to $LOGFILE"