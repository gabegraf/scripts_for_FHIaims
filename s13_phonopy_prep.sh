#!/bin/bash

for file in geometry.in-*; do 
	number=$(echo "$file" | grep -oP '(?<=geometry.in-)\d+')

	folder="${number}-displacement"

	mkdir -p "$folder"
	mv "$file" "$folder/geometry.in"

	cp jobscript.sh control.in "$folder/"

	cd "$folder" && sbatch jobscript.sh && cd ../
done

echo "All jobs have been submitted!"
