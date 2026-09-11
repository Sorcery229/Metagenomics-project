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
reproducibility/  Backing data for the manuscript's low-similarity analysis
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

## Running the notebooks

Notebook file paths are relative to the repository root, so launch Jupyter
from the repo root:

```
pip install -r requirements.txt
jupyter notebook
```

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

All references that point to files committed in this repository have been
updated to their new locations under `data/`, `results/`, and `figures/`.

## Reproducibility: low-similarity cases

`reproducibility/` backs the manuscript analysis of cross-database pairs below
50% similarity:

- `low_similarity_461_classification.csv` — every low-similarity pair with its
  accessions, similarity value, source-repository `genome_status`, and an
  assigned category. Summing the `category` column reproduces the reported
  counts: 461 pairs total, of which 7 were already flagged by the source
  repository (5 deprecated + 2 with no status) and 199 arise from a
  plasmid-vs-complete retrieval-convention difference.
- `low_similarity_461_accessions.csv` — the accession list for those pairs.

See `reproducibility/DATA_DICTIONARY.md` for column definitions.
