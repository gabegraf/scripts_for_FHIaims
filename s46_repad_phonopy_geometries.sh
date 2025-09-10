#!/bin/bash

for f in geometry.in-*; do
    # Extract the number after the dash
    num=${f#geometry.in-}

    # Force decimal (avoids octal problems)
    num=$((10#$num))

    # Pad to 4 digits
    newnum=$(printf "%04d" "$num")

    # New filename
    newname="geometry.in-${newnum}"

    # Rename if different
    if [ "$f" != "$newname" ]; then
        mv "$f" "$newname"
        echo "Renamed $f -> $newname"
    fi
done

