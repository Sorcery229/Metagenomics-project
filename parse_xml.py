import csv
import pickle
import time  # To generate unique timestamps
from xml.etree import ElementTree as ET
from collections import Counter  # Import Counter to count occurrences

# Hardcoded subject length for the single subject (can be overwritten by CSV file)
def get_subject_length(subject_length_from_csv=None):
    if subject_length_from_csv:
        return subject_length_from_csv  # Use subject length from the CSV file if provided
    return 1600000  # Default hardcoded length for the subject

# Step 2: Parse the BLAST XML file to extract hits and HSPs
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

                for hsp in hit.findall(".//Hsp"):
                    hsp_data = {
                        "Hsp_qseq": hsp.find("Hsp_qseq").text,
                        "Hsp_hseq": hsp.find("Hsp_hseq").text,
                        "Hsp_hit-from": hsp.find("Hsp_hit-from").text,
                        "Hsp_hit-to": hsp.find("Hsp_hit-to").text,
                        "Hsp_align-len": hsp.find("Hsp_align-len").text,
                    }
                    hit_data["Hit_hsps"].append(hsp_data)

                hits.append(hit_data)
        print(f"Parsed {len(hits)} hits.")
    except Exception as e:
        print(f"Error parsing XML: {e}")
    
    return hits

# Step 3: Create the subject array
def create_subject_array(subject_length, hit_list):
    # Initialize the subject array with zeros for one big subject
    subject_array = [0] * subject_length

    # Process each hit and HSP
    for hit in hit_list:
        for hsp in hit["Hit_hsps"]:
            qseq = hsp["Hsp_qseq"]
            hseq = hsp["Hsp_hseq"]
            hit_from = int(hsp["Hsp_hit-from"])
            hit_to = int(hsp["Hsp_hit-to"])
            alignment_len = int(hsp["Hsp_align-len"])

            # Determine if it's a reverse match based on hit_from and hit_to
            reverse = hit_from > hit_to
            start, end = (hit_to, hit_from) if reverse else (hit_from, hit_to)
            print(start,end)
            
            for i in range(alignment_len):
                #print(alignment_len)
                if reverse:
                    seq_pos = end - 1 - i
                    q_char = qseq[alignment_len - i - 1]
                    h_char = hseq[alignment_len - i - 1]
                else:
                    seq_pos = start - 2 + i
                    q_char = qseq[i]
                    h_char = hseq[i]

                # Match conditions (1 and 5 for normal/reverse match, 2 and 6 for gaps)
                if q_char == h_char and q_char != '-':
                    #print(seq_pos)

                    subject_array[seq_pos] = 1 if not reverse else 5  # Normal match → 1, Reverse match → 5
                elif h_char == '-':
                    subject_array[seq_pos] = 2 if not reverse else 6  # Gap in hseq → 2, Reverse gap in hseq → 6
                elif q_char == '-':
                    subject_array[seq_pos] = 3 if not reverse else 7  # Gap in qseq → 3, Reverse gap in qseq → 7
                else:
                    subject_array[seq_pos] = 4  # Mismatch → 4 for both normal and reverse

    print("Subject array created.")
    return subject_array

# Step 4: Save the subject array to a binary file using pickle
def save_subject_array_to_binary(subject_array, output_filename):
    try:
        with open(output_filename, 'wb') as file:
            pickle.dump(subject_array, file)  # Save array as binary using pickle
        print(f"Subject array saved to binary file: {output_filename}")
    except Exception as e:
        print(f"Error saving subject array to binary: {e}")

# Step 5: Count occurrences of each value (n of 0, n of 1, etc.)
def count_value_occurrences(subject_array):
    counts = Counter(subject_array)  # Count occurrences of each number
    count_dict = {}
    for number in range(8):  # Assuming values 0 to 7
        count_dict[f'n_of_{number}'] = counts.get(number, 0)
    return count_dict

# Function to process the CSV file
def process_csv(input_csv, output_csv):
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        # Define fieldnames explicitly, ensuring all new columns are included
        fieldnames = reader.fieldnames + ['n_of_0', 'n_of_1', 'n_of_2', 'n_of_3', 'n_of_4', 'n_of_5', 'n_of_6', 'n_of_7', 'binary_filename']  # New columns for counts and binary filename
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            # Take subject length from the "subject_length" column and XML filename from the "blast_output" column
            subject_length = int(row['subject_length'])  # "subject_length" column
            xml_file = row['blast_output']  # "blast_output" column

            hit_list = parse_hits_from_xml(xml_file)

            if hit_list:
                subject_array = create_subject_array(subject_length, hit_list)
                value_counts = count_value_occurrences(subject_array)

                # Add the counts to the row
                for number, count in value_counts.items():
                    row[number] = count

                # Generate a unique binary filename using the current timestamp
                timestamp = int(time.time())
                binary_filename = f"{row['blast_output'].split('.')[0]}_subject_array_{timestamp}.bin"
                save_subject_array_to_binary(subject_array, binary_filename)

                # Add the binary filename to the row
                row['binary_filename'] = binary_filename

                writer.writerow(row)  # Write row with new counts and binary filename
            else:
                row.update({f'n_of_{i}': 0 for i in range(8)})  # Set all counts to 0
                row['binary_filename'] = ''  # No binary file if no hits
                writer.writerow(row)

# Example usage
input_csv = 'genomes.csv'  # Replace with your input CSV file
output_csv = 'output_file.csv'  # Replace with desired output CSV file

process_csv(input_csv, output_csv)

print("Processing complete. Counts saved to output_file.csv and subject arrays saved to binary files.")
