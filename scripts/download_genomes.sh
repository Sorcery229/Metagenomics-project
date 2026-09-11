#!/bin/bash

# Input CSV file
datafile="genomes.csv"

# Count total number of genomes to process
total_genomes=$(($(wc -l < "$datafile") - 1))
processed_genomes=0

echo "Downloading genomes... Total: $total_genomes"

# Read the CSV file line by line (excluding the header)
tail -n +2 "$datafile" | while IFS=, read -r line; do
    # Extract genome_id and ftp_path correctly while ignoring other columns
    genome_id=$(echo "$line" | cut -d, -f4)
    ftp_path=$(echo "$line" | cut -d, -f6)
    
    processed_genomes=$((processed_genomes + 1))
    remaining_genomes=$((total_genomes - processed_genomes))
    echo "Processing genome: $genome_id ($processed_genomes/$total_genomes, $remaining_genomes left)"
    
    # Construct the first FTP download URL
    bvbrc_url="ftp://ftp.bvbrc.org/genomes/$genome_id/$genome_id.fna"
    bvbrc_filename="$genome_id.fna"
    
    # Download the file from BVBRC
    echo "Downloading from BVBRC: $bvbrc_url"
    curl -# "$bvbrc_url" -o "$bvbrc_filename"
    
    # Construct the second FTP download URL
    ncbi_filename="$(basename "$ftp_path")_genomic.fna.gz"
    ncbi_url="$ftp_path/$ncbi_filename"
    
    # Download and unzip the file from NCBI
    echo "Downloading from NCBI: $ncbi_url"
    curl -# "$ncbi_url" -o "$ncbi_filename"
    gunzip -f "$ncbi_filename"
    ncbi_unzipped_filename="${ncbi_filename%.gz}"
    
    echo "Downloaded and extracted: $ncbi_unzipped_filename"
    
    # Append filenames to the original CSV by modifying the same file
    # Use `awk` to add the filenames to the original CSV in the correct place
    echo "$line,$bvbrc_filename,$ncbi_unzipped_filename" >> temp_file.csv

done

# Replace the original CSV with the updated one
mv temp_file.csv "$datafile"

echo "Download complete. Original file updated."

