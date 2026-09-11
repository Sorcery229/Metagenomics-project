#!/bin/bash

# Get the directory where the script is located
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SUMMARY="$DIR/summary_with_filenames.tsv"
FNA_SOURCE="$DIR/virushostdb.genomic.fna"
TEMP="$DIR/updated_summary.tsv"
BLAST="$HOME/Desktop/gboldirev1/blast/bin/blastn"  # <-- updated BLAST path

# Header for the updated summary
echo -e "Index\tSequence ID\tVirus Tax ID\tFTP Path\tRefSeq Filename\tQuery Length\tSubject Length\tBLAST Output" > "$TEMP"

# Function to extract a single indexed sequence from a multi-FASTA file
extract_fasta_by_index() {
    local index="$1"
    local input="$2"
    local output="$3"
    awk -v idx="$index" '
        BEGIN {found=0; count=0}
        /^>/ {
            count++
            if (count == idx) {
                print; found=1; next
            }
            if (count > idx && found == 1) {
                exit
            }
        }
        found == 1 {print}
    ' "$input" > "$output"
}

# Read the summary file line by line (skip header)
tail -n +2 "$SUMMARY" | while IFS=$'\t' read -r index seqid taxid ftp refseq_filename; do
    # Remove carriage returns or other whitespace
    refseq_filename=$(echo "$refseq_filename" | tr -d '\r')
    seqid=$(echo "$seqid" | tr -d '\r')

    echo "Processing index $index ($seqid)..."

    # Full path to RefSeq file
    refseq_path="$DIR/$refseq_filename"

    # Check that the subject (RefSeq) file exists
    if [[ ! -f "$refseq_path" ]]; then
        echo "Missing RefSeq file: $refseq_path"
        continue
    fi

    # Extract query fasta (from .fna source)
    query_fasta="$DIR/query_${index}.fna"
    extract_fasta_by_index "$index" "$FNA_SOURCE" "$query_fasta"

    # Ensure file is not empty
    if [[ ! -s "$query_fasta" ]]; then
        echo "Empty or missing query sequence at index $index"
        continue
    fi

    # Flatten sequences but preserve header
    for f in "$query_fasta" "$refseq_path"; do
        header=$(grep "^>" "$f")
        sequence=$(grep -v "^>" "$f" | tr -d '\n')
        echo -e "$header\n$sequence" > "$f"
    done

    # Get lengths
    query_length=$(grep -v "^>" "$query_fasta" | tr -d '\n' | wc -c)
    subject_length=$(grep -v "^>" "$refseq_path" | tr -d '\n' | wc -c)

    # Define output file
    blast_output="$DIR/${index}_${seqid}_blast.xml"

    # Run BLAST
    "$BLAST" -query "$query_fasta" -subject "$refseq_path" -out "$blast_output" -outfmt 5 \
        || echo "BLAST failed for $query_fasta vs $refseq_path" >> "$DIR/failed_blast_jobs.log"

    # Write output row
    echo -e "$index\t$seqid\t$taxid\t$ftp\t$refseq_filename\t$query_length\t$subject_length\t$blast_output" >> "$TEMP"
done

echo "BLAST complete. Summary written to $TEMP"
