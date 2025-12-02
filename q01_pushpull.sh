#!/usr/bin/env bash

# This is a personal script for organizing my filesystem on my home cluster
# and syncing with remote HPC resources. It is not intended for public use.

set -euo pipefail
IFS=$'\n\t'

# Usage: q01_cluster_pushpull.sh [push|pull] [projects|cluster] [kestrel|dcc|perlmutter] [folder]
if [ "$#" -eq 4 ]; then
    action="$1"
    location="$2"
    computer="$3"
    folder="$4"
elif [ "$#" -eq 3 ]; then
    action="$1"
    location="$2"
    computer="$3"
    folder=""
else
    echo "Usage: $0 [push|pull] [projects|cluster] [kestrel|dcc|perlmutter] [folder]"
    exit 1
fi

case "$action" in
    push|pull) ;;
    *) echo "Action must be 'push' or 'pull'"; exit 1 ;;
esac

case "$location" in
    cluster|projects) ;;
    *) echo "Location must be 'cluster' or 'projects'"; exit 1 ;;
esac

: "${username:?username must be set in the environment}"
: "${timewarp_id:?timewarp_id must be set in the environment}"
: "${kestrel_id:?kestrel_id must be set in the environment}"
: "${dcc_id:?dcc_id must be set in the environment}"
: "${perlmutter_id:?perlmutter_id must be set in the environment}"

case "$computer" in
    kestrel)
        computer_id="$kestrel_id"
        active_dir="02_kestrel"
        remote_dir="/$username"
        ;;
    dcc)
        computer_id="$dcc_id"
        active_dir="01_dcc"
        remote_dir="/work/$username"
        ;;
    perlmutter)
        computer_id="$perlmutter_id"
        active_dir="03_perlmutter"
        remote_dir="/pscratch/sd/g/$username/03_perlmutter"
        ;;
    *) echo "Unknown computer: $computer"; exit 1 ;;
esac

command -v rsync >/dev/null 2>&1 || { echo "rsync not found"; exit 1; }
command -v globus >/dev/null 2>&1 || echo "Warning: globus not found; cluster transfers will fail"

folder="${folder%/}"

if [ "$location" = "cluster" ]; then
    if [ "$action" = "push" ]; then
        if [ -z "$folder" ]; then
            src="$timewarp_id:/home/$username/active/$active_dir"
            dest="$computer_id:$remote_dir/$active_dir"
        else
            src="$timewarp_id:/home/$username/active/$active_dir/$folder"
            dest="$computer_id:$remote_dir/$folder"
        fi
        globus transfer --recursive --sync-level checksum "$src" "$dest"
    else
        # pull
        globus transfer --sync-level checksum --recursive --exclude "*.csc" \
            --delete-destination-extra --skip-source-errors \
            "$computer_id:$remote_dir/" "$timewarp_id:/home/$username/active/$active_dir/"
    fi

elif [ "$location" = "projects" ]; then
    if [ -z "$folder" ]; then
        echo "When using projects location you must specify a folder"
        exit 1
    fi

    project_src="$HOME/projects/$folder"
    basefolder="$(basename "$folder")"
    active_dest="$HOME/active/$active_dir/$basefolder"

    if [ "$action" = "pull" ]; then
        if [ ! -d "$project_src" ]; then
            echo "Source project does not exist: $project_src"; exit 1
        fi
        rsync -a "$project_src/" "$active_dest/"
        readlink -f "$project_src" > "$active_dest/timewarp_path.txt"
    else
        # push (move back to original location)
        tp="$HOME/active/$active_dir/$folder/timewarp_path.txt"
        if [ ! -f "$tp" ]; then
            echo "Missing timewarp_path.txt in $HOME/active/$active_dir/$folder"; exit 1
        fi
        src_abs="$(< "$tp")"
        if [ -z "$src_abs" ]; then
            echo "timewarp_path.txt is empty"; exit 1
        fi

        # ensure destination exists
        mkdir -p "$src_abs"

        rsync -a --delete --remove-source-files --exclude 'timewarp_path.txt' \
            "$HOME/active/$active_dir/$folder/" "$src_abs/"

        rm -f "$tp"
        find "$HOME/active/$active_dir/$folder" -depth -type d -empty -exec rmdir {} + || true
    fi
fi




