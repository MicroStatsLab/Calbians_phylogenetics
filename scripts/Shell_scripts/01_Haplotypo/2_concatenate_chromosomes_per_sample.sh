#!/bin/bash

# Loop through all files with a .fasta extension
for file in ../*.fasta; do
    # Extract the base name before the first underscore
    base_name=$(basename "$file" | cut -d'.' -f1)
    
    # Define the output file name based on the base name
    output_file="${base_name}.fasta"

    echo "Processing $file -> $output_file..."

    # Concatenate sequences within the file, ignoring lines starting with ">"
    awk '/^>/ {next} {printf "%s", $0} END {print ""}' "$file" > "$output_file"
    
    # Add the new header based on the base name
    sed -i "1i >${base_name}" "$output_file"
done

echo "All files processed."
