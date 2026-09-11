#!/bin/bash

datafile="accessions_with_fasta.csv"
tempfile="temp_accessions_with_fasta.csv"

# Extract header and find the index of "ftp_path"
header=$(head -n 1 "$datafile")
IFS=',' read -r -a columns <<< "$header"

ftp_index=-1
for i in "${!columns[@]}"; do
    if [[ "${columns[$i]}" == "ftp_path" ]]; then
        ftp_index=$((i + 1))  # CSV fields are 1-based for `cut`
        break
    fi
done

if [[ $ftp_index -eq -1 ]]; then
    echo "Error: ftp_path column not found in header."
    exit 1
fi

# Write updated header
echo "$header,downloaded_filename" > "$tempfile"

# Process each line (skip header)
tail -n +2 "$datafile" | while IFS= read -r line; do
    # Extract ftp_path using the correct column index
    ftp_path=$(echo "$line" | cut -d',' -f"$ftp_index")

    # Skip empty ftp_path
    if [[ -z "$ftp_path" ]]; then
        echo "$line," >> "$tempfile"
        continue
    fi

    # Build filename and URL
    filename="$(basename "$ftp_path")_genomic.fna.gz"
    url="$ftp_path/$filename"

    echo "Downloading: $url"

    # Download file
    curl -# "$url" -o "$filename"

    # Append filename to line
    echo "$line,$filename" >> "$tempfile"
done

# Replace original file with updated file
mv "$tempfile" "$datafile"

echo "Download complete. Updated file: $datafile"
