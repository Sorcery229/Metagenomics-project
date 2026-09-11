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
                        # If subject (hseq) has a gap, encode insertion in query at the last valid subject position
                        if prev_valid_pos is not None and subject_array[prev_valid_pos] > 0:
                            subject_array[prev_valid_pos] = -subject_array[prev_valid_pos]  # Make it negative to indicate insertion
                        continue  # Skip modifying this position

                    if 0 <= seq_pos < subject_length:
                        if q_char == h_char and q_char != '-':
                            subject_array[seq_pos] = 1 if not reverse else 5
                        elif q_char == '-':
                            subject_array[seq_pos] = 3 if not reverse else 7  # Insertion in subject
                        else:
                            subject_array[seq_pos] = 4  # Mismatch

                        filled_positions.add(seq_pos)
                        prev_valid_pos = seq_pos  # Update last valid subject position

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
