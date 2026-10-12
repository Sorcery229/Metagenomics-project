# Worked example / smoke test

A minimal, self-contained test of the BLAST-XML parsing stage
(`scripts/parse_xml.py`), the step that turns a BLAST alignment into the
per-base similarity counts used throughout the study. It needs only Python
(standard library) - no BLAST install and no network.

## Inputs
- `example_blast.xml` - a single BLAST hit: a 20 bp subject aligned to a query
  that is identical except for two mismatches (at positions 5 and 15).
- `genomes.csv` - the driver table `parse_xml.py` reads. Required columns:
  `subject_length` and `blast_output` (path to the XML); any other columns are
  carried through to the output.

## Run
```
cd example
python ../scripts/parse_xml.py
```
This writes `output_file.csv` (and a pickled `*_subject_array_*.bin`).

## Expected output
Compare the count columns of `output_file.csv` against `expected_output.csv`:

| column  | value | meaning |
|---------|-------|---------|
| n_of_1  | 18    | matched bases |
| n_of_4  | 2     | mismatched bases |
| all other n_of_* | 0 | no gaps/insertions in this example |

Per-base similarity = (n_of_1 + n_of_5) / subject_length = 18 / 20 = 0.90.

`expected_output.csv` omits the `binary_filename` column, whose value contains
a run timestamp and therefore differs on every run; only the `n_of_*` count
columns are deterministic and should be compared.

## One-line check
```
cd example && python ../scripts/parse_xml.py \
  && python -c "import csv; r=list(csv.DictReader(open('output_file.csv')))[0]; \
     assert (r['n_of_1'],r['n_of_4'])==('18','2'), r; print('OK: similarity =', (int(r['n_of_1'])+int(r['n_of_5']))/int(r['subject_length']))"
```


## Figure reproduction

`reproduce_similarity_figure.py` regenerates the cross-database per-base
similarity distribution from committed data (`results/output_file.csv`), using
the same binning as `notebooks/visualization.ipynb`:

```
python example/reproduce_similarity_figure.py
```

Outputs `similarity_distribution.png` and `similarity_counts.csv`. Diff the
counts against the committed `expected_similarity_counts.csv` (20,651 pairs;
17,633 at 100%, 2,238 at 99%, 70 below 10%, 33 below 70%). Counts are
deterministic; the PNG may vary at the pixel level across matplotlib versions.
