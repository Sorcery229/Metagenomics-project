# Cross-DB Genomic Comparator (CDGC)

A framework for quantifying genome-similarity discrepancies across microbial
reference databases (BV-BRC/PATRIC, RefSeq, Ensembl). Genome pairs matched
across databases are aligned and compared at base level (matches, mismatches,
insertions, deletions) to produce a per-pair similarity value, which is then
used to flag inconsistencies between the databases.

## Repository layout

```
notebooks/        Analysis and figure notebooks
  compare_genomes.ipynb     pairwise comparison of matched genomes
  match_databases.ipynb     cross-database strain matching
  stats.ipynb               summary statistics
  visualization.ipynb       similarity-distribution figures
scripts/          Pipeline scripts
  download_genomes.sh       fetch genome FASTAs from BV-BRC and NCBI
  run_blast.sh              run pairwise BLAST over a genome list
  parse_xml.py              parse BLAST XML into per-pair similarity CSV
data/             Input reference data
  assembly_summaries/       NCBI assembly_summary tables (fungi, viral)
  genomes_fungi.csv         fungal genome list
  updated_summary.tsv       consolidated summary table
  contigs/                  multi-contig genome file lists
results/          Per-pair similarity output tables (output_file*.csv, 2to2results.csv)
figures/          Generated figures (output_viral.png)
databases/        Self-contained per-database pipelines (each with its own
                  download/blast scripts, genome list, and outputs)
  refseq_bvbrc/  refseq_ensembl_fungi_data/  refseq_fungi_db/  refseq_viral_host_db/
example/          Worked examples that run against committed data:
  parse_xml smoke test and a figure-reproduction script
```

## Pipeline overview

1. `scripts/download_genomes.sh` reads a genome list (`genomes.csv`) and
   downloads the matching FASTA from BV-BRC and NCBI.
2. `scripts/run_blast.sh` runs pairwise BLAST between the paired genomes.
3. `scripts/parse_xml.py` parses the BLAST XML output into a per-pair
   similarity table (`output_file.csv`).
4. The `notebooks/` read those result tables to match strains across
   databases, compute statistics, and render the similarity-distribution
   figures.

The four folders under `databases/` are stand-alone reruns of this same
pipeline against individual database pairs; each uses only paths local to its
own folder and can be run from inside that folder.

## Installation

The recommended way to reproduce the full environment, including the aligner,
is the committed conda environment (`environment.yml`), which installs the
pinned Python stack together with NCBI BLAST+ 2.11.0 and puts `blastn` on your
`PATH`:

```
git clone https://github.com/Sorcery229/Metagenomics-project.git
cd Metagenomics-project
conda env create -f environment.yml
conda activate cdgc
```

Alternatively, install the Python dependencies with pip and supply BLAST+
yourself:

```
pip install -r requirements.txt
```

The alignment step (`scripts/run_blast.sh`) requires the NCBI BLAST+ `blastn`
executable, which is not pip-installable. The analysis in this study used
**NCBI BLAST+ 2.11.0**, which the conda environment above pins. If you install
BLAST+ manually (from https://blast.ncbi.nlm.nih.gov/), confirm the version with
`blastn -version`.

### Software versions

Verified with **Python 3.11**: pandas 2.3.3, numpy 2.4.6, matplotlib 3.11.0,
seaborn 0.13.2, matplotlib-venn 1.1.1 (pinned in `requirements.txt`), and
**NCBI BLAST+ 2.11.0** (pinned in `environment.yml`).

### Input and output formats

- `scripts/download_genomes.sh` - input: a genome list `genomes.csv`
  (one genome id per row); output: downloaded FASTA (`.fna`) files.
- `scripts/run_blast.sh` - input: paired FASTA files; output: one BLAST XML
  report per pair.
- `scripts/parse_xml.py` - input: a CSV with `subject_length` and
  `blast_output` (path to a BLAST XML) columns; output: `output_file.csv`, the
  same rows augmented with per-base count columns `n_of_-7` .. `n_of_7`
  (1 = match, 4 = mismatch, 3/7 = insertion in subject, negative = insertion in
  query, 0 = uncovered) plus a `binary_filename` column. Per-base similarity of
  a pair is `(n_of_1 + n_of_5) / subject_length`.

## Running the pipeline (exact commands)

Each script reads and writes fixed file names in the current working directory,
so the "input parameters" are the files that must be present in that directory,
not command-line flags. Run each step from a directory that contains the named
input.

1. **Download genomes.** Requires `genomes.csv` in the working directory, with
   the BV-BRC genome id in column 4 and the NCBI FTP path in column 6:

   ```
   bash scripts/download_genomes.sh
   ```
   Output: one `<genome_id>.fna` (BV-BRC) and one `*_genomic.fna` (NCBI,
   gunzipped) per row.

2. **Pairwise alignment.** Requires `genomes.csv` with the paired FASTA paths in
   columns 8 (query) and 9 (subject). Either activate the conda environment (so
   `blastn` is on `PATH`) or edit the `BLAST` variable at the top of the script
   to point at your `blastn`:

   ```
   bash scripts/run_blast.sh
   ```
   Each pair is aligned with the full-length command:
   ```
   blastn -query <query.fna> -subject <subject.fna> -out <pair>.xml -outfmt 5
   ```
   Output: one BLAST XML report (`-outfmt 5`) per pair, and `genomes.csv`
   updated with `query_length`, `subject_length`, and `blast_output` columns.

3. **Parse alignments.** Requires `genomes.csv` with `subject_length` and
   `blast_output` columns (as produced by step 2) in the working directory:

   ```
   python scripts/parse_xml.py
   ```
   Reads `genomes.csv`, writes `output_file.csv` (the per-base count columns
   `n_of_-7` .. `n_of_7` plus `binary_filename`). Per-pair similarity =
   `(n_of_1 + n_of_5) / subject_length`.

The four per-database pipelines under `databases/` repeat these same three steps
with the same fixed file names (`genomes.csv` -> `output_file.csv`); run each
from inside its own folder. Each `databases/*/run_blast*.sh` has its own
hard-coded `BLAST` path at the top that must be edited to your `blastn` before
running.

The figure-reproduction scripts take no inputs beyond the committed data and are
run from the repository root:

```
python example/reproduce_similarity_figure.py    # -> example/similarity_distribution.png, example/similarity_counts.csv
python scripts/reproduce_similarity_figures.py   # -> figures/fig_similarity_{strain,species}.png, figures/similarity_band_distribution.csv
python scripts/reproduce_venn_overlaps.py        # -> figures/venn_{strain,species}_overlap.png, figures/venn_overlap_counts.csv
```

## Worked example (smoke test)

`example/` contains a minimal, network-free test of the parsing stage that
runs with only the Python standard library:

```
cd example
python ../scripts/parse_xml.py
```

Compare the `n_of_*` columns of the generated `output_file.csv` against the
committed `example/expected_output.csv` (`n_of_1 = 18`, `n_of_4 = 2`, giving a
similarity of 18/20 = 0.90). See `example/README.md` for the one-line check.

## Reproducing the notebooks

Notebook file paths are relative to the repository root, so launch Jupyter
from the repo root:

```
jupyter notebook
```

`notebooks/visualization.ipynb` renders the similarity-distribution figures
(for example `figures/output_viral.png`, the viral-genome similarity
distribution); `notebooks/compare_genomes.ipynb` and `notebooks/stats.ipynb`
produce the pairwise comparisons and summary statistics; `notebooks/match_databases.ipynb`
performs the cross-database strain matching. Note that some notebook cells
depend on the external inputs listed below.

### External inputs not tracked in this repo

The notebooks were originally authored on the first author's workstation and
still reference several inputs that are not committed here. These must be
supplied locally before the affected cells will run:

- `genomes.csv` — the working genome list consumed/produced by the scripts
- a bacterial `assembly_summary.txt` (NCBI RefSeq bacteria summary)
- `species_EnsemblBacteria.txt`, `refseq_representative_assemblies.csv`,
  `patric_representative_assemblies.csv` — representative-assembly lists
- `summary.csv`, `10k_diff.csv` and the intermediate files generated during a
  run (`remaining_genomes.csv`, `common_strains_with_*.csv`,
  `complete_labeled_genomes.csv`, `all_bacteria_results.csv`)

For full reproducibility we recommend depositing these external inputs (the
working genome lists and the large representative-assembly tables that are too
big to track in git) in a Zenodo record and citing its DOI here, so that a
reader can retrieve the exact files used rather than re-deriving them. Several
of these inputs come from living databases, so the download date or release of
each file matters; record it alongside the deposited files (see the acquisition
dates noted for the committed snapshots in the Venn section below).

All references that point to files committed in this repository have been
updated to their new locations under `data/`, `results/`, and `figures/`.

## Reproducing a main figure

`example/reproduce_similarity_figure.py` regenerates the cross-database
per-base similarity distribution directly from committed data
(`results/output_file.csv`), using the same binning as
`notebooks/visualization.ipynb` (per-pair similarity = `n_of_1 / subject_length`):

```
python example/reproduce_similarity_figure.py
```

It writes `example/similarity_distribution.png` (the bar chart) and
`example/similarity_counts.csv` (category -> count). Diff the counts against
the committed `example/expected_similarity_counts.csv` to confirm the result:
20,651 pairs, 17,633 at 100% identity and 2,238 at 99%, with 70 below 10% and
33 below 70%. The counts are deterministic; the PNG may differ at the pixel
level across matplotlib versions.

## Cross-database per-base similarity figures

`scripts/reproduce_similarity_figures.py` regenerates the two similarity-distribution
figures from the committed per-domain tables under `data/similarity/`:

```
python scripts/reproduce_similarity_figures.py
```

Inputs (`data/similarity/`), one row per cross-database genome pair:

| column          | description                                                        |
|-----------------|--------------------------------------------------------------------|
| `subject_length`| length of the reference (subject) genome, in bases                  |
| `matched_bases` | aligned identical bases (bacteria: `n_of_1 + n_of_5`; viruses/fungi: aligned matches) |
| `similarity`    | `matched_bases / subject_length`                                    |
| `source`        | (fungi only) `ensembl` or `fungidb`                                 |

The script asserts `similarity == matched_bases / subject_length` on load, so the
values are recomputed from the primitive counts rather than trusted as stored.

Similarity bands (identical for all three domains):

- `100%`     : `similarity >= 0.999999`
- `95-<100%` : `0.95 <= similarity < 0.999999`
- `<95%`     : `similarity < 0.95`

Outputs (`figures/`):

- `fig_similarity_strain.png`  — Bacteria, strain-level (N = 285,051)
- `fig_similarity_species.png` — Viruses (N = 12,715) and Fungi (N = 470), species-level
- `similarity_band_distribution.csv` — the per-band counts, percentages, mean and
  standard deviation behind both figures

Expected band percentages (also in `figures/similarity_band_distribution.csv`):

| domain   | level  | N       | 100%  | 95-<100% | <95%  | mean   | sd     |
|----------|--------|---------|-------|----------|-------|--------|--------|
| Bacteria | strain | 285,051 | 55.7% | 42.1%    | 2.3%  | 0.9947 | 0.0443 |
| Viruses  | species| 12,715  | 98.7% | 0.7%     | 0.6%  | 0.9972 | 0.0411 |
| Fungi    | species| 470     | 15.3% | 74.7%    | 10.0% | 0.9807 | 0.0608 |

The counts are deterministic; the PNGs may differ at the pixel level across
matplotlib versions.

## Database overlap Venn figures (strain and species level)

`scripts/reproduce_venn_overlaps.py` regenerates two Venn figures directly from
committed database metadata; every region count is a set operation, nothing is
hard-coded:

```
python scripts/reproduce_venn_overlaps.py
```

It writes `figures/venn_strain_overlap.png` (Bacteria: RefSeq vs BV-BRC),
`figures/venn_species_overlap.png` (Viruses and Fungi), and
`figures/venn_overlap_counts.csv` (region -> reproduced count alongside the
manuscript value).

Set definitions (key -> committed source file):

| Domain | Level | Database | Key | Source file |
| --- | --- | --- | --- | --- |
| Bacteria | strain | RefSeq | (taxid, strain) | `data/venn/bacteria_refseq_strain_keys.csv.gz` |
| Bacteria | strain | BV-BRC | (taxon_id, strain) | `data/venn/bacteria_bvbrc_strain_keys.csv.gz` |
| Viruses | species | RefSeq | taxid | `data/assembly_summaries/assembly_summary_viral.txt` |
| Viruses | species | VirusHostDB | virus tax id | `data/venn/virushostdb.tsv.gz` |
| Fungi | species | Ensembl | taxonomy_id | `data/venn/species_EnsemblFungi.txt` |
| Fungi | species | RefSeq | taxid | `data/assembly_summaries/assembly_summary_fungi_refseq.txt` |
| Fungi | species | FungiDB | species_ncbi_tax_id | `data/venn/fungidb_organisms.csv` |

The bacterial strain-key files are the deduplicated (taxid, strain) representative
keys derived from NCBI `assembly_summary.txt` (RefSeq representative genomes) and
the BV-BRC `genome_metadata` table; the raw tables are too large to commit.
`virushostdb.tsv.gz` and `fungidb_organisms.csv` are snapshots of the living
VirusHostDB (https://www.genome.jp/ftp/db/virushostdb/) and FungiDB
(https://fungidb.org) lists. Because those two databases keep growing, the
VirusHostDB-only and FungiDB region counts drift slightly above the manuscript
values, while the RefSeq/BV-BRC strain intersection and the Ensembl/RefSeq fungal
intersection reproduce the manuscript exactly.

**Acquisition dates of the committed snapshots** (record these when citing the
data, since these are living databases): the VirusHostDB list
(https://www.genome.jp/ftp/db/virushostdb/) and the FungiDB organism list
(https://fungidb.org, VEuPathDB AllOrganisms report) were downloaded on
2026-09-15; the Ensembl Fungi species list
(`data/venn/species_EnsemblFungi.txt`) was obtained on 2026-09-15; and the NCBI
RefSeq assembly summaries (`data/assembly_summaries/assembly_summary_viral.txt`
and `assembly_summary_fungi_refseq.txt`) were downloaded on 2026-09-11.
