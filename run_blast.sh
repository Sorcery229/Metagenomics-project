#!/bin/bash

# Set input CSV file
CSV_FILE="genomes.csv"
TEMP_CSV="temp_file.csv"

# Path to BLAST tool
BLAST="/Users/gboldirev1/blast/bin/blastn"

# Create a temporary CSV file and add the header
echo "col1,col2,col3,col4,col5,col6,col7,query_fna,subject_fna,query_length,subject_length,blast_output" > "$TEMP_CSV"

# Function to process and concatenate sequences in an .fna file
process_fna() {
    local fna_file="$1"
    local seq_length

    # Check if the file contains multiple sequences using grep ">"
    if grep -q "^>" "$fna_file"; then
        # Concatenate all sequences into a single sequence
        awk '/^>/ {if (seq) print seq; seq=""} !/^>/ {seq=seq$0} END {if (seq) print seq}' "$fna_file" > "${fna_file}.tmp"

        # Replace the original file with the processed one
        mv "${fna_file}.tmp" "$fna_file"
    fi

    # Calculate sequence length after concatenation
    seq_length=$(grep -v "^>" "$fna_file" | tr -d '\n' | wc -c)
    echo "$seq_length"
}

# Get the total number of lines in the CSV file
total_lines=$(wc -l < "$CSV_FILE")
current_line=0

# Read CSV line by line and process each row
while IFS=, read -r col1 col2 col3 col4 col5 col6 col7 query_fna subject_fna rest; do
    ((current_line++))
    
    # Display progress
    echo "Processing $current_line of $total_lines..."

    if [[ -f "$query_fna" && -f "$subject_fna" ]]; then
        # Process and concatenate sequences in both files
        query_length=$(process_fna "$query_fna")
        subject_length=$(process_fna "$subject_fna")

        # Generate a unique BLAST output filename
        blast_output="${query_fna//\//_}_vs_${subject_fna//\//_}_blast.xml"

        # Run BLAST
        "$BLAST" -query "$query_fna" -subject "$subject_fna" -out "$blast_output" -outfmt 5

        # Overwrite the row in the temp file
        echo "$col1,$col2,$col3,$col4,$col5,$col6,$col7,$query_fna,$subject_fna,$query_length,$subject_length,$blast_output" >> "$TEMP_CSV"
    else
        # Copy the original row if files are missing
        echo "$col1,$col2,$col3,$col4,$col5,$col6,$col7,$query_fna,$subject_fna,$rest" >> "$TEMP_CSV"
    fi
done < "$CSV_FILE"

# Replace the original CSV with the updated one
mv "$TEMP_CSV" "$CSV_FILE"
echo "Processing complete. CSV file updated."
