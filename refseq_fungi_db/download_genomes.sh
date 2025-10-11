#!/bin/bash
set -euo pipefail

# Always execute in the folder where the job was submitted
cd "${SGE_O_WORKDIR:-$(pwd)}"

# -------- Config --------
datafile="genomes.csv"
download_dir="downloads"
temp_csv="temp_file.csv"
header="index,taxid,strain_name,genome_id,contig_count_fungidb,ftp_path,fasta_link,query_fna,subject_fna"

# Check that download directory exists
if [[ ! -d "$download_dir" ]]; then
  echo "ERROR: download directory '$download_dir' does not exist."
  exit 1
fi

: > download_errors.log

# -------- Downloader (no curl) --------
have_wget=0
command -v wget >/dev/null 2>&1 && have_wget=1

download_file() {
  # usage: download_file URL OUT_PATH
  local url="$1"
  local out="$2"
  [[ -z "$url" ]] && return 1

  if [[ $have_wget -eq 1 ]]; then
    wget -q --tries=3 --timeout=60 -O "$out" "$url"
    return $?
  else
    # Python stdlib fallback, with retries
    python3 - "$url" "$out" <<'PY'
import sys, time, urllib.request
if len(sys.argv) < 3:
    sys.exit(1)
url, out = sys.argv[1], sys.argv[2]
for attempt in range(3):
    try:
        req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as r, open(out, "wb") as f:
            f.write(r.read())
        sys.exit(0)
    except Exception as e:
        if attempt < 2:
            time.sleep(2 * (attempt+1))
        else:
            sys.stderr.write(str(e) + "\n")
            sys.exit(1)
PY
  fi
}

# -------- CSV writer helper (proper quoting) --------
csv_escape() {
  local s=${1//$'\r'/}
  s=${s//$'\n'/ }
  s=${s//\"/\"\"}
  printf '"%s"' "$s"
}
append_csv_row() {
  {
    csv_escape "$1"; printf ','
    csv_escape "$2"; printf ','
    csv_escape "$3"; printf ','
    csv_escape "$4"; printf ','
    csv_escape "$5"; printf ','
    csv_escape "$6"; printf ','
    csv_escape "$7"; printf ','
    csv_escape "$8"; printf ','
    csv_escape "$9"; printf '\n'
  } >> "$temp_csv"
}

# -------- Count & progress --------
if [[ ! -f "$datafile" ]]; then
  echo "ERROR: $datafile not found." >&2
  exit 1
fi

total_genomes=$(( $(wc -l < "$datafile") - 1 ))
processed_genomes=0
echo "Downloading genomes... Total: $total_genomes"

# Start temp CSV with header
echo "$header" > "$temp_csv"

# -------- Main loop --------
while IFS=$'\t' read -r \
  index taxid strain_name genome_id contig_count_fungidb ftp_path fasta_link
do
  [[ -z "${genome_id// }" ]] && continue

  processed_genomes=$((processed_genomes + 1))
  remaining_genomes=$((total_genomes - processed_genomes))
  echo "Processing genome: $genome_id ($processed_genomes/$total_genomes, $remaining_genomes left)"

  ftp_path=${ftp_path//\"/}
  fasta_link=${fasta_link//\"/}

  [[ -z "$fasta_link" ]] && echo "Warning: empty FungiDB link for $genome_id" >> download_errors.log
  [[ -z "$ftp_path"   ]] && echo "Warning: empty NCBI path for $genome_id"   >> download_errors.log

  # --- FungiDB download ---
  query_fna=""
  if [[ -n "$fasta_link" ]]; then
    fasta_filename="$(basename "$fasta_link")"
    fasta_out="$download_dir/$fasta_filename"
    echo "Downloading FungiDB: $fasta_link"
    if download_file "$fasta_link" "$fasta_out"; then
      query_fna="$fasta_out"
    else
      echo "Failed: $fasta_link" >> download_errors.log
    fi
  fi

  # --- NCBI / RefSeq download ---
  subject_fna=""
  if [[ -n "$ftp_path" ]]; then
    ftp_path="${ftp_path%/}"
    base="$(basename "$ftp_path")"
    ncbi_url="$ftp_path/${base}_genomic.fna.gz"
    ncbi_filename="${base}_genomic.fna.gz"
    ncbi_out="$download_dir/$ncbi_filename"

    echo "Downloading NCBI: $ncbi_url"
    if download_file "$ncbi_url" "$ncbi_out"; then
      if [[ -s "$ncbi_out" ]]; then
        if gunzip -f "$ncbi_out"; then
          subject_fna="$download_dir/${base}_genomic.fna"
          echo "Downloaded and extracted: $subject_fna"
        else
          echo "gunzip failed (kept gz): $ncbi_out" >> download_errors.log
          subject_fna="$ncbi_out"
        fi
      else
        echo "Empty or missing after download: $ncbi_out" >> download_errors.log
      fi
    else
      echo "Failed: $ncbi_url" >> download_errors.log
    fi
  fi

  append_csv_row "$index" "$taxid" "$strain_name" "$genome_id" "$contig_count_fungidb" "$ftp_path" "$fasta_link" "$query_fna" "$subject_fna"

done < <(python3 - "$datafile" <<'PY'
import sys, csv
path = sys.argv[1]
with open(path, newline='') as f:
    r = csv.DictReader(f)
    cols = ["index","taxid","strain_name","genome_id","contig_count_fungidb","ftp_path","fasta_link"]
    for row in r:
        out = [row.get(c, "") for c in cols]
        print("\t".join(out))
PY
)

mv "$temp_csv" "$datafile"

echo "Download complete. Files saved in '$download_dir'."
echo "Check download_errors.log for any failed downloads."

