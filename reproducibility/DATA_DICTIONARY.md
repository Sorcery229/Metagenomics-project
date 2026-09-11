# Data dictionary — low-similarity classification

`low_similarity_461_classification.csv` — one row per cross-database genome
pair with a per-base similarity below 50% (461 rows).

| column                | description |
|-----------------------|-------------|
| genome_id             | BV-BRC/PATRIC genome identifier used to join the pair |
| name                  | organism name |
| refseq_acc            | RefSeq-side accession |
| bv_assembly_accession | BV-BRC-side assembly accession |
| subject_length        | length of the RefSeq-side sequence |
| bv_len                | length of the BV-BRC-side sequence |
| sim                   | per-base similarity value for the pair (0-1) |
| same_assembly         | whether both records map to the same underlying assembly |
| bv_genome_status      | source-repository status from the BV-BRC `genome_status` field |
| category              | assigned class (see below) |

## Categories (sum of `category` reproduces every reported count)

| category | count | meaning |
|----------|-------|---------|
| valid genome - discrepancy surfaced by CDGC | 230 | genuine cross-database sequence difference newly surfaced by CDGC |
| plasmid (retrieval convention) | 199 | a plasmid record from one database compared against the complete sequence from the other |
| no current BV-BRC record | 25 | genome_id no longer served by the source repository |
| flagged (deprecated) | 5 | source repository marks the assembly deprecated |
| flagged (missing status) | 2 | source repository record has no assigned status |

Summary used in the manuscript and reviewer reply: **461 total** = **7 already
flagged by the source** (5 deprecated + 2 missing status) + **199 plasmid
retrieval-convention** + **254 not previously flagged** (230 genuine differences
surfaced by CDGC + 25 with no current source record).

`low_similarity_461_accessions.csv` — the accession list for the same 461 pairs.
