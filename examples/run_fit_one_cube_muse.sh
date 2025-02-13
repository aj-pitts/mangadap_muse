#!/bin/bash

ARGS="$@"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

LOGDIR="DAP_logs"
LOGSUBDIR="$LOGDIR/$1_$2"
LOGFILE="$LOGSUBDIR/DAP_$1_$2_${TIMESTAMP}.log"

mkdir -p "$LOGSUBDIR"

nohup python fit_one_cube_muse.py "$@" > "$LOGFILE" 2>&1 &

echo "Running Script: python fit_one_cube_muse.py "$@" > "$LOGFILE" 2>&1 &"
echo "Logging output to $LOGFILE"