# SeqCarver

Carve any genome into feature and non-feature sequence blocks from a FASTA + GFF3/GTF annotation. Outputs a FASTA and a coordinate table covering 100% of the input sequence.

---

## What it does

SeqCarver walks every sequence in your genome file from position 1 to the end and slices it into two kinds of blocks:

- **Feature blocks** — regions that match a feature type in your annotation (default: `gene`). Overlapping features are merged into a single block.
- **Non-feature blocks** — everything in between, labeled dynamically as `non_<feature-type>` (e.g. `non_gene`, `non_exon`, `non_CDS`).

Every single base in the input ends up in exactly one block — nothing is dropped.

## Outputs

| File | Contents |
|------|----------|
| `<stem>_<feature-type>_SeqCarver.fasta` | One FASTA record per block, in genome order |
| `<stem>_<feature-type>_SeqCarver.tsv` | Coordinates, block type, length, strand, name, biotype, locus tag |

## Requirements

- Python 3.7+
- [pyfaidx](https://github.com/mdshw5/pyfaidx)

```bash
pip install pyfaidx
```

## Usage

```bash
# basic -- output written to current directory
python3 SeqCarver.py --seq genome.fasta --annotation genes.gff3

# write to a specific output directory
python3 SeqCarver.py --seq genome.fasta --annotation genes.gff3 --out-dir results/

# carve by exon instead of gene
python3 SeqCarver.py --seq genome.fasta --annotation genes.gff3 --feature-type exon

# gzipped inputs, gzipped FASTA output
python3 SeqCarver.py --seq genome.fa.gz --annotation genes.gtf.gz --gzip-fasta
```

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--seq` | required | Sequence file: `.fasta` / `.fa` / `.fna` / `.fastq` / `.fq` (`.gz` OK) |
| `--annotation` | required | Annotation file: `.gff` / `.gff3` / `.gtf` (`.gz` OK) |
| `--feature-type` | `gene` | GFF/GTF column-3 type to use as feature blocks |
| `--out-dir` | `.` | Directory to write output files into |
| `--out-prefix` | auto | Override auto-naming with a full path prefix |
| `--gzip-fasta` | off | Gzip the output FASTA |
| `--workdir` | auto temp | Scratch dir for decompressed/converted files |

## Notes

- If no annotation rows match `--feature-type`, the script lists every feature type it did find so you can pick the right one.
- Sequences with no annotated features produce a single non-feature block spanning the whole sequence.
- Both GFF3 and GTF formats are auto-detected; `.gz` compression is handled transparently.
