import csv
import pickle
import time  # To generate unique timestamps
from xml.etree import ElementTree as ET
from collections import Counter  # Import Counter to count occurrences

def parse_hits_from_xml(xml_file):
    hits = []
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        for iteration in root.findall(".//Iteration"):
            iteration_hits = iteration.find("Iteration_hits")
            for hit in iteration_hits.findall("Hit"):
                hit_data = {
                    "Hit_num": hit.find("Hit_num").text,
                    "Hit_id": hit.find("Hit_id").text,
                    "Hit_len": hit.find("Hit_len").text,
                    "Hit_hsps": []
                }

                hsps = []
                for hsp in hit.findall(".//Hsp"):
                    hsp_data = {
                        "Hsp_qseq": hsp.find("Hsp_qseq").text,
                        "Hsp_hseq": hsp.find("Hsp_hseq").text,
                        "Hsp_hit-from": int(hsp.find("Hsp_hit-from").text) - 1,  # Convert to 0-based
                        "Hsp_hit-to": int(hsp.find("Hsp_hit-to").text) - 1,  # Convert to 0-based
                        "Hsp_align-len": int(hsp.find("Hsp_align-len").text),
                        "Hsp_identity": float(hsp.find("Hsp_identity").text)  # Extract identity
                    }
                    hsps.append(hsp_data)

                hit_data["Hit_hsps"] = sorted(hsps, key=lambda x: (-x["Hsp_identity"], -x["Hsp_align-len"]))
                hits.append(hit_data)
        print(f"Parsed {len(hits)} hits from {xml_file}.")
    except Exception as e:
        print(f"Error parsing XML ({xml_file}): {e}")
    
    return hits

def create_subject_array(subject_length, hit_list):
    subject_array = [0] * subject_length  # Initialize subject array
    filled_positions = set()

    for hit in hit_list:
        for hsp in hit["Hit_hsps"]:
            try:
                qseq = hsp["Hsp_qseq"]
                hseq = hsp["Hsp_hseq"]
                hit_from = hsp["Hsp_hit-from"]
                hit_to = hsp["Hsp_hit-to"]
                alignment_len = hsp["Hsp_align-len"]

                reverse = hit_from > hit_to
                start, end = (hit_to, hit_from) if reverse else (hit_from, hit_to)
                prev_valid_pos = None  # Keep track of last valid subject position

                for i in range(alignment_len):
                    if reverse:
                        seq_pos = end - i
                        if seq_pos in filled_positions or not (0 <= seq_pos < subject_length):
                            continue
                        q_char = qseq[alignment_len - i - 1]
                        h_char = hseq[alignment_len - i - 1]
                    else:
                        seq_pos = start + i
                        if seq_pos in filled_positions or not (0 <= seq_pos < subject_length):
                            continue
                        q_char = qseq[i]
                        h_char = hseq[i]

                    if h_char == '-':  
                        if prev_valid_pos is not None and subject_array[prev_valid_pos] > 0:
                            subject_array[prev_valid_pos] = -subject_array[prev_valid_pos]  # Mark insertion
                        continue

                    if 0 <= seq_pos < subject_length:
                        if q_char == h_char and q_char != '-':
                            subject_array[seq_pos] = 1 if not reverse else 5
                        elif q_char == '-':
                            subject_array[seq_pos] = 3 if not reverse else 7  # Insertion in subject
                        else:
                            subject_array[seq_pos] = 4  # Mismatch

                        filled_positions.add(seq_pos)
                        prev_valid_pos = seq_pos

            except Exception as e:
                print(f"Error processing HSP: {e}")

    print("Subject array created.")
    return subject_array

def save_subject_array_to_binary(subject_array, output_filename):
    try:
        with open(output_filename, 'wb') as file:
            pickle.dump(subject_array, file)
        print(f"Subject array saved to: {output_filename}")
    except Exception as e:
        print(f"Error saving subject array to binary: {e}")

def count_value_occurrences(subject_array):
    counts = Counter(subject_array)
    return {f'n_of_{num}': counts.get(num, 0) for num in range(-7, 8)}  # Include negative values

def process_csv(input_csv, output_csv):
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + [f'n_of_{i}' for i in range(-7, 8)] + ['binary_filename']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            try:
                subject_length_str = row.get('subject_length')
                if not subject_length_str or subject_length_str.strip() == '':
                    print(f"Skipping row due to missing subject_length: {row}")
                    continue

                subject_length = int(float(subject_length_str.strip().replace(',', '')))
                xml_file = row['blast_output']
                hit_list = parse_hits_from_xml(xml_file)

                if hit_list:
                    subject_array = create_subject_array(subject_length, hit_list)
                    value_counts = count_value_occurrences(subject_array)

                    for number, count in value_counts.items():
                        row[number] = count

                    timestamp = int(time.time())
                    binary_filename = f"{row['blast_output'].split('.')[0]}_subject_array_{timestamp}.bin"
                    save_subject_array_to_binary(subject_array, binary_filename)

                    row['binary_filename'] = binary_filename
                    writer.writerow(row)
                else:
                    print(f"No hits found in {xml_file}, skipping.")

            except Exception as e:
                print(f"Skipping row due to error: {e}")

input_csv = 'genomes.csv'
output_csv = 'output_file.csv'

process_csv(input_csv, output_csv)

print("Processing complete. Valid results saved to output_file.csv, binary files saved.")

[akarlsbe@login2 databases] cat run_blast.sh
#!/bin/bash

# safer, but we handle errors manually
set -u
cd /u/home/a/akarlsbe/scratch/databases || exit 1

# absolute path to BLAST
BLAST="/u/local/apps/blast+/2.11.0/ncbi-blast-2.11.0+/bin/blastn"

# log helper
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

# setup files
CSV_FILE="genomes.csv"
TEMP_CSV="temp_file.csv"
LOG_FILE="blast_debug.log"

echo "===== Starting BLAST job =====" > "$LOG_FILE"
log "Working directory: $(pwd)"

# check BLAST
if [[ ! -x "$BLAST" ]]; then
    log "ERROR: blastn not found at $BLAST. Exiting."
    exit 1
fi
log "Using BLAST binary at $BLAST"

# check CSV
if [[ ! -f "$CSV_FILE" ]]; then
    log "ERROR: CSV file not found: $CSV_FILE"
    exit 1
fi

# prepare CSV header
echo "col1,col2,col3,col4,col5,col6,col7,query_fna,subject_fna,query_length,subject_length,blast_output" > "$TEMP_CSV"

total_lines=$(($(wc -l < "$CSV_FILE") - 1))
log "Total comparisons to process: $total_lines"

current_line=0

# main loop
tail -n +2 "$CSV_FILE" | while IFS=, read -r index taxid strain_name genome_id contig_count_fungidb ftp_path fasta_link query_fna subject_fna; do
    ((current_line++))
    log "---- [${current_line}/${total_lines}] ----"
    log "Query: $query_fna | Subject: $subject_fna"

    # cleanup quotes
    query_fna=$(echo "$query_fna" | tr -d '"')
    subject_fna=$(echo "$subject_fna" | tr -d '"')

    # check files
    if [[ ! -f "$query_fna" ]]; then
        log "WARNING: Query file missing: $query_fna"
        continue
    fi
    if [[ ! -f "$subject_fna" ]]; then
        log "WARNING: Subject file missing: $subject_fna"
        continue
    fi

    # compute lengths
    query_len=$(grep -v "^>" "$query_fna" | tr -d '\n' | wc -c || echo 0)
    subject_len=$(grep -v "^>" "$subject_fna" | tr -d '\n' | wc -c || echo 0)
    log "Query length: $query_len | Subject length: $subject_len"

    # output filename
    output_xml="$(basename "${query_fna%.*}")_vs_$(basename "${subject_fna%.*}")_blast.xml"

    # run blast with error handling
    log "Running: $BLAST -query $query_fna -subject $subject_fna -out $output_xml -outfmt 5"
    if "$BLAST" -query "$query_fna" -subject "$subject_fna" -out "$output_xml" -outfmt 5 2>>"$LOG_FILE"; then
        log "✅ BLAST success for $query_fna vs $subject_fna"
    else
        code=$?
        log "❌ BLAST failed (exit code $code) for $query_fna vs $subject_fna"
        continue
    fi

    # record result
    echo "$index,$taxid,$strain_name,$genome_id,$contig_count_fungidb,$ftp_path,$fasta_link,$query_fna,$subject_fna,$query_len,$subject_len,$output_xml" >> "$TEMP_CSV"

done

log "Merging temporary file into $CSV_FILE"
mv "$TEMP_CSV" "$CSV_FILE"
log "All done — BLAST processing finished successfully."

