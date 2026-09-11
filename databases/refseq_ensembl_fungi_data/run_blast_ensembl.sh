#!/bin/bash

CSV_FILE="genomes.csv"
TEMP_CSV="temp_file.csv"
BLAST="/Users/gboldirev1/Desktop/gboldirev1/blast/bin/blastn"

# Write header including the empty first column
echo ",division,assembly,bioproject,organism_name,infraspecific_name,#assembly_accession,genome_size,fasta_filename,downloaded_filename,query_length,subject_length,blast_output" > "$TEMP_CSV"

# Function to process fasta files and return sequence length
process_fna() {
    local fna_file="$1"
    echo "Processing fasta file: $fna_file" >&2

    if grep -q "^>" "$fna_file"; then
        echo "Concatenating sequences by removing headers..." >&2
        sed '/^>/d' "$fna_file" > "${fna_file}.tmp"
        mv "${fna_file}.tmp" "$fna_file"
    else
        echo "No headers found, skipping concatenation." >&2
    fi

    local seq_length
    seq_length=$(grep -v "^>" "$fna_file" | tr -d '\n' | wc -c)
    echo "$seq_length"
}

total_lines=$(wc -l < "$CSV_FILE")
echo "Total lines in CSV (including header): $total_lines"

current_line=0
exec 3< "$CSV_FILE"

# Read and skip header
read -r header <&3
echo "CSV header: $header"

while IFS=',' read -r unused division assembly bioproject organism_name infraspecific_name assembly_accession genome_size fasta_filename downloaded_filename <&3; do
    ((current_line++))
    echo "----------------------------------" >&2
    echo "Processing row $current_line of $((total_lines - 1))" >&2
    echo "unused='$unused'" >&2
    echo "division='$division'" >&2
    echo "assembly='$assembly'" >&2
    echo "bioproject='$bioproject'" >&2
    echo "organism_name='$organism_name'" >&2
    echo "infraspecific_name='$infraspecific_name'" >&2
    echo "assembly_accession='$assembly_accession'" >&2
    echo "genome_size='$genome_size'" >&2
    echo "fasta_filename='$fasta_filename'" >&2
    echo "downloaded_filename='$downloaded_filename'" >&2

    fasta_filename=$(echo "$fasta_filename" | xargs)
    downloaded_filename=$(echo "$downloaded_filename" | xargs)

    if [[ -z "$fasta_filename" || -z "$downloaded_filename" ]]; then
        echo "Empty fasta or downloaded filename, skipping row." >&2
        continue
    fi

    if [[ "$downloaded_filename" == *.gz ]]; then
        subject_unzipped="${downloaded_filename%.gz}"
        echo "Unzipping $downloaded_filename to $subject_unzipped" >&2
        if gunzip -c "$downloaded_filename" > "$subject_unzipped"; then
            echo "Unzip successful." >&2
        else
            echo "Failed to unzip $downloaded_filename. Skipping row." >&2
            continue
        fi
    else
        subject_unzipped="$downloaded_filename"
        echo "Subject file is not gzipped. Using as is." >&2
    fi

    if [[ ! -f "$fasta_filename" ]]; then
        echo "Query fasta file missing: $fasta_filename. Skipping row." >&2
        continue
    fi

    if [[ ! -f "$subject_unzipped" ]]; then
        echo "Subject fasta file missing: $subject_unzipped. Skipping row." >&2
        continue
    fi

    query_length=$(process_fna "$fasta_filename")
    subject_length=$(process_fna "$subject_unzipped")

    blast_output="${fasta_filename//\//_}_vs_${subject_unzipped//\//_}_blast.xml"
    echo "Running BLAST: query=$fasta_filename subject=$subject_unzipped" >&2
    if "$BLAST" -query "$fasta_filename" -subject "$subject_unzipped" -out "$blast_output" -outfmt 5; then
        echo "BLAST finished successfully." >&2
    else
        echo "BLAST failed." >&2
    fi

    echo "Writing results to output CSV" >&2
    echo ",$division,$assembly,$bioproject,$organism_name,$infraspecific_name,$assembly_accession,$genome_size,$fasta_filename,$downloaded_filename,$query_length,$subject_length,$blast_output" >> "$TEMP_CSV"

done

exec 3<&-

mv "$TEMP_CSV" "$CSV_FILE"
echo "All rows processed. Updated CSV saved to $CSV_FILE." >&2
