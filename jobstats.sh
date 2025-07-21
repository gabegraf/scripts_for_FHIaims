#!/usr/bin/env bash

# A script containing functions to parse Slurm job statistics.
# Provides a useful summary of seff, sstat, and sacct.
#
# Usage: jobstats.sh <job_id> [n_steps_for_sacct (optional)]
#
# Authors: Uthpala Herath and ChatGPT

# A function to parse Slurm strings like "7721840K", "5.73M", "186.30M", etc.
# into a byte count. Also handles decimal numbers plus an optional K/M/G suffix.
to_bytes() {
    local val="$1"
    # If empty or not recognized, return 0.
    [[ -z "$val" ]] && echo 0 && return

    # Regex: optionally a decimal number, then optionally K/M/G
    # Examples that match:
    #   "12345" -> no suffix, assume raw bytes
    #   "123.45M"
    #   "186.30M"
    #   "11521316K"
    #   "5.73M"
    if [[ "$val" =~ ^([0-9]+(\.[0-9]+)?)([KMG])?$ ]]; then
        local num="${BASH_REMATCH[1]}"   # e.g. "186.30"
        local unit="${BASH_REMATCH[3]}"  # e.g. "M"
        local factor=1

        case "$unit" in
            K) factor=$((1024)) ;;
            M) factor=$((1024*1024)) ;;
            G) factor=$((1024*1024*1024)) ;;
            *) factor=1 ;;  # no suffix => bytes
        esac

        # We must do a floating-point multiply: num * factor.
        # Use awk to produce an integer result (rounding).
        awk -v n="$num" -v f="$factor" 'BEGIN {
            bytes = n * f
            # round to nearest integer
            printf "%.0f", bytes
        }'
    else
        # If purely integer digits, assume raw bytes
        if [[ "$val" =~ ^[0-9]+$ ]]; then
            echo "$val"
        else
            # Unrecognized format
            echo 0
        fi
    fi
}

# Convert bytes to human-readable MB or GB (switch at 1GB).
to_mb_or_gb() {
    local bytes="$1"
    local oneGB=$((1024*1024*1024))
    if (( bytes < oneGB )); then
        # Show in MB
        awk -v b="$bytes" 'BEGIN { printf "%.2fMB", b/(1024*1024) }'
    else
        # Show in GB
        awk -v b="$bytes" 'BEGIN { printf "%.2fGB", b/(1024*1024*1024) }'
    fi
}

# The main function:
#    - Parse "Nodes: X" from seff output
#    - seff summary
#    - sstat summary (JobID,JobName,MaxRSS,MaxDiskWrite) in MB/GB, plus "Total MaxRSS"
#    - sacct summary (JobID,JobName,MaxRSS,MaxDiskWrite) in MB/GB, plus "Total MaxRSS"
jobstats() {
    local jobid="$1"
    local maxSteps="$2"   # optional second argument

    if [[ -z "$jobid" ]]; then
        echo "Usage: jobstats.sh <job_id> [<n_steps>]"
        return 1
    fi

    # (A) seff output
    local seff_output
    seff_output="$(seff "$jobid" 2>&1)"
    echo "=== [1/3] seff summary ==="
    echo "$seff_output"
    echo

    # Parse "Nodes: X" from seff, if present
    local numNodes
    numNodes="$(echo "$seff_output" | sed -n 's/^Nodes:\s*\([0-9]\+\).*/\1/p')"
    [[ -z "$numNodes" ]] && numNodes=1

    # (B) sstat summary
    echo "=== [2/3] sstat summary (LIVE) ==="
    local sstat_out
    sstat_out="$(
        sstat --noheader --format=JobID%30,MaxRSS,MaxDiskWrite \
              -j "${jobid}" 2>/dev/null
    )"

    if [[ -z "$sstat_out" ]]; then
        echo "Job ${jobid} is not running. Check sacct summary."
    else
        echo "JobID step      | MaxRSS/node | Total MaxRSS | MaxDiskWrite"
        echo "----------------|------------ | ------------ | ------------"
        while IFS= read -r line; do
            [[ -z "$line" ]] && continue
            IFS=' ' read -r stepJobID rawRss rawDisk <<< "$line"
            [[ -z "$stepJobID" ]] && continue

            local rssBytes diskBytes
            rssBytes="$(to_bytes "$rawRss")"
            diskBytes="$(to_bytes "$rawDisk")"

            local rssHuman diskHuman
            rssHuman="$(to_mb_or_gb "$rssBytes")"
            diskHuman="$(to_mb_or_gb "$diskBytes")"

            local totalRSSBytes=$(( rssBytes * numNodes ))
            local totalRSSHuman
            totalRSSHuman="$(to_mb_or_gb "$totalRSSBytes")"

            printf "%-15s | %-11s | %-12s | %s\n" \
                "$stepJobID" "$rssHuman" "$totalRSSHuman" "$diskHuman"
        done <<< "$sstat_out"
    fi
    echo

    # (C) sacct summary
    echo "=== [3/3] sacct summary (WARNING: Inactive until job step completion!) ==="
    local sacct_out
    sacct_out="$(
        sacct --noheader \
              --format=JobID%30,JobName%25,MaxRSS,MaxDiskWrite \
              -j "${jobid}" 2>/dev/null
    )"

    if [[ -z "$sacct_out" ]]; then
        echo "No sacct info found."
    else
        # Convert sacct_out into an array so we can optionally truncate
        mapfile -t sacct_lines <<< "$sacct_out"

        # If maxSteps is given and numeric, we tail only last N lines
        if [[ -n "$maxSteps" && "$maxSteps" =~ ^[0-9]+$ ]]; then
            if (( maxSteps < ${#sacct_lines[@]} )); then
                # Truncate to last N lines
                sacct_lines=("${sacct_lines[@]: -$maxSteps}")
            fi
        fi

        echo "JobID step      | JobName                    | MaxRSS/node | Total MaxRSS | MaxDiskWrite"
        echo "----------------|----------------------------|-------------|--------------|-------------"
        for line in "${sacct_lines[@]}"; do
            [[ -z "$line" ]] && continue
            IFS=' ' read -r cJobID cJobName cMaxRSS cMaxDiskWrite <<< "$line"
            [[ -z "$cJobID" ]] && continue

            local rssBytes diskBytes
            rssBytes="$(to_bytes "$cMaxRSS")"
            diskBytes="$(to_bytes "$cMaxDiskWrite")"

            local rssHuman diskHuman
            rssHuman="$(to_mb_or_gb "$rssBytes")"
            diskHuman="$(to_mb_or_gb "$diskBytes")"

            local totalRSSBytes=$(( rssBytes * numNodes ))
            local totalRSSHuman
            totalRSSHuman="$(to_mb_or_gb "$totalRSSBytes")"

            printf "%-15s | %-26s | %-11s | %-12s | %s\n" \
                "$cJobID" "$cJobName" "$rssHuman" "$totalRSSHuman" "$diskHuman"
        done
    fi
}

# Call the function with any arguments passed to the script
jobstats "$@"

