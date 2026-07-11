#!/usr/bin/env python3
"""
SeqCarver -- carve a genome into feature and non-feature blocks.

Walks every sequence from position 1 to end and splits it into:
  - feature blocks     : annotation rows matching --feature-type (overlaps merged)
  - non-feature blocks : everything in between (labeled non_<feature-type>,
                         e.g. non_gene, non_exon, non_CDS)

The label is dynamic on purpose -- "intergenic" is only meaningful when
--feature-type is "gene". Run with --feature-type exon and those gaps are
non_exon, which is what they actually are.

Supported inputs:
  sequences   : .fasta / .fa / .fna / .fastq / .fq  (plain or .gz)
  annotations : .gff / .gff3 / .gtf                 (plain or .gz)

Every base in the input ends up in exactly one output block -- nothing dropped.

Outputs:
  <out-dir>/<stem>_<feature-type>_SeqCarver.fasta   one record per block, genome order
  <out-dir>/<stem>_<feature-type>_SeqCarver.tsv     coords / type / length / name / biotype / strand

Usage:
  python3 SeqCarver.py --seq genome.fasta --annotation genes.gff3
  python3 SeqCarver.py --seq genome.fasta --annotation genes.gff3 --out-dir results
  python3 SeqCarver.py --seq genome.fasta --annotation genes.gff3 --feature-type exon
  python3 SeqCarver.py --seq genome.fa.gz --annotation genes.gtf --gzip-fasta

If no rows match --feature-type the script prints what feature types it did
find and exits, so you can pick the right one.
"""

import argparse
import csv
import gzip
import os
import re
import shutil
import sys
import tempfile

try:
    from pyfaidx import Fasta
except ImportError:
    sys.exit("Missing dependency. Install with: pip install pyfaidx --break-system-packages")


# attribute regexes -- GFF3 uses key=value, GTF uses key "value"
GFF3_ATTR_RE = re.compile(r'([\w\.]+)=([^;]*)')
GTF_ATTR_RE  = re.compile(r'([\w\.]+)\s+"([^"]*)"')


def sniff_annotation_format(path, opener):
    """Peek at the first non-comment, non-blank feature line to decide gff3 vs gtf."""
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = cols[8]
            if GTF_ATTR_RE.search(attrs) and "=" not in attrs.split(";")[0]:
                return "gtf"
            if GFF3_ATTR_RE.search(attrs):
                return "gff3"
            return "gtf"  # fallback
    return "gff3"


def parse_attrs(attr_str, fmt):
    if fmt == "gtf":
        return dict(GTF_ATTR_RE.findall(attr_str))
    return dict(GFF3_ATTR_RE.findall(attr_str))


# attribute key priority lists -- covers GFF3/GTF/Ensembl/NCBI/GENCODE
NAME_KEYS    = ["Name", "gene_name", "gene", "transcript_name", "ID", "gene_id"]
ID_KEYS      = ["ID", "gene_id", "transcript_id"]
BIOTYPE_KEYS = ["gene_biotype", "biotype", "gene_type", "transcript_biotype", "transcript_type"]
LOCUS_KEYS   = ["locus_tag", "old_locus_tag", "Alias"]


def first_present(d, keys, default=""):
    for k in keys:
        if k in d and d[k]:
            return d[k]
    return default



def opener_for(path):
    return gzip.open if path.endswith(".gz") else open


def is_fastq(path):
    p = path[:-3] if path.endswith(".gz") else path
    return p.lower().endswith((".fastq", ".fq"))


def ensure_fasta(seq_path, workdir):
    """Get an uncompressed FASTA in workdir that pyfaidx can index.

    Handles .gz decompression and FASTQ -> FASTA conversion as needed.
    pyfaidx writes a .fai next to the file, so it needs to live somewhere writable.
    """
    base = os.path.basename(seq_path)
    if base.endswith(".gz"):
        base = base[:-3]

    if is_fastq(seq_path):
        out_path = os.path.join(workdir, os.path.splitext(base)[0] + ".fasta")
        print(f"Converting FASTQ -> FASTA: {out_path}")
        opener = opener_for(seq_path)
        with opener(seq_path, "rt") as fin, open(out_path, "w") as fout:
            while True:
                header = fin.readline()
                if not header:
                    break
                seq = fin.readline()
                plus = fin.readline()
                qual = fin.readline()
                if not seq:
                    break
                name = header.strip().lstrip("@")
                fout.write(f">{name}\n{seq.strip()}\n")
        return out_path

    # plain or gzipped FASTA -- copy to workdir so the .fai can be written there
    out_path = os.path.join(workdir, base)
    opener = opener_for(seq_path)
    if seq_path.endswith(".gz"):
        print(f"Decompressing {seq_path} -> {out_path}")
        with opener(seq_path, "rt") as fin, open(out_path, "w") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copy(seq_path, out_path)
    return out_path


def derive_stem(seq_path):
    """Pull a clean base name out of a sequence file path.

    /data/GCF_000001215.fna.gz -> GCF_000001215
    genome.fastq               -> genome
    """
    name = os.path.basename(seq_path)
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in (".fasta", ".fa", ".fna", ".fastq", ".fq"):
        if name.lower().endswith(ext):
            name = name[:-len(ext)]
            break
    return name or "output"



def load_features_by_seq(annotation_path, feature_type):
    opener = opener_for(annotation_path)
    fmt = sniff_annotation_format(annotation_path, opener)
    print(f"Annotation format detected: {fmt.upper()}")

    features_by_seq = {}
    seen_types = {}
    with opener(annotation_path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols[:9]
            seen_types[ftype] = seen_types.get(ftype, 0) + 1
            if ftype.lower() != feature_type.lower():
                continue
            a = parse_attrs(attrs, fmt)
            name = first_present(a, NAME_KEYS, default="unknown")
            feat_id = first_present(a, ID_KEYS, default="")
            biotype = first_present(a, BIOTYPE_KEYS, default="NA")
            locus_tag = first_present(a, LOCUS_KEYS, default="")
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            if start_i > end_i:
                start_i, end_i = end_i, start_i
            features_by_seq.setdefault(seqid, []).append(
                (start_i, end_i, strand, name, feat_id, biotype, locus_tag)
            )

    if not features_by_seq:
        print(f"\nWARNING: no rows with feature type '{feature_type}' found.")
        print("Feature types present in this annotation file:")
        for t, c in sorted(seen_types.items(), key=lambda x: -x[1]):
            print(f"    {t}: {c}")
        sys.exit(1)

    for seqid in features_by_seq:
        features_by_seq[seqid].sort(key=lambda x: (x[0], x[1]))
    return features_by_seq


def merge_overlaps(feature_list):
    """Merge overlapping/touching intervals, keeping all member features in each block."""
    merged = []
    for feat in feature_list:
        start, end, strand, name, feat_id, biotype, locus_tag = feat
        if merged and start <= merged[-1]["end"]:
            merged[-1]["end"] = max(merged[-1]["end"], end)
            merged[-1]["features"].append(feat)
        else:
            merged.append({"start": start, "end": end, "features": [feat]})
    return merged



def write_fasta_record(fh, header, seq_str, width=70):
    fh.write(">" + header + "\n")
    for i in range(0, len(seq_str), width):
        fh.write(seq_str[i:i + width] + "\n")


def main():
    print("=== SeqCarver ===")
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--seq",          required=True, help="Sequence file: .fasta/.fa/.fna/.fastq/.fq (.gz OK)")
    ap.add_argument("--annotation",   required=True, help="Annotation file: .gff/.gff3/.gtf (.gz OK)")
    ap.add_argument("--out-dir",      default=".",   help="Output directory (default: current dir)")
    ap.add_argument("--out-prefix",   default=None,
                    help="Override auto-naming with a full path prefix "
                         "(e.g. results/run1 -> results/run1.fasta + .tsv)")
    ap.add_argument("--feature-type", default="gene",
                    help="GFF/GTF column-3 type to treat as feature blocks (default: gene)")
    ap.add_argument("--gzip-fasta",   action="store_true", help="Gzip the output FASTA")
    ap.add_argument("--workdir",      default=None,
                    help="Scratch dir for decompressed/converted files (default: auto temp dir)")
    args = ap.parse_args()

    workdir = args.workdir or tempfile.mkdtemp(prefix="geneextract_")
    os.makedirs(workdir, exist_ok=True)

    if args.out_prefix is None:
        stem = derive_stem(args.seq)
        filename = f"{stem}_{args.feature_type}_SeqCarver"
        args.out_prefix = os.path.join(args.out_dir, filename)
        print(f"--out-prefix not given, using: {args.out_prefix}")

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print("Preparing sequence file...")
    fasta_path = ensure_fasta(args.seq, workdir)

    print("Loading annotation...")
    features_by_seq = load_features_by_seq(args.annotation, args.feature_type)
    n_feat = sum(len(v) for v in features_by_seq.values())
    print(f"  {n_feat} '{args.feature_type}' features across {len(features_by_seq)} sequences")

    print("Indexing sequence file...")
    fasta = Fasta(fasta_path, sequence_always_upper=False)

    out_fasta_path = args.out_prefix + ".fasta"
    out_table_path = args.out_prefix + ".tsv"

    fasta_out = open(out_fasta_path, "w")
    table_out = open(out_table_path, "w", newline="")
    tsv = csv.writer(table_out, delimiter="\t")
    tsv.writerow(["seqid", "block_type", "start", "end", "length", "strand",
                  "name", "feature_id", "biotype", "locus_tag", "n_features_merged"])

    n_seqs = len(fasta.keys())
    total_feat_blocks = 0
    total_other_blocks = 0
    non_feature_label = f"non_{args.feature_type}"

    for i, seqid in enumerate(fasta.keys(), 1):
        seq_len = len(fasta[seqid])
        feat_list = features_by_seq.get(seqid, [])
        blocks = merge_overlaps(feat_list) if feat_list else []

        # walk the sequence left-to-right, interleaving feature and non-feature segments
        cursor   = 1
        segments = []
        for b in blocks:
            if b["start"] > cursor:
                segments.append((non_feature_label, cursor, b["start"] - 1, None))
            segments.append((args.feature_type, b["start"], b["end"], b))
            cursor = b["end"] + 1
        if cursor <= seq_len:
            segments.append((non_feature_label, cursor, seq_len, None))
        if not segments and seq_len > 0:
            segments = [(non_feature_label, 1, seq_len, None)]

        for seg_type, start, end, block in segments:
            length = end - start + 1
            if length <= 0:
                continue
            seq_str = str(fasta[seqid][start - 1:end])

            if seg_type == non_feature_label:
                header = f"{seqid}:{start}-{end}|{non_feature_label}"
                write_fasta_record(fasta_out, header, seq_str)
                tsv.writerow([seqid, non_feature_label, start, end, length, ".", "", "", "", "", 0])
                total_other_blocks += 1
            else:
                members = block["features"]
                names = ";".join(m[3] for m in members)
                ids = ";".join(m[4] for m in members if m[4])
                biotypes = ";".join(sorted(set(m[5] for m in members)))
                locus_tags = ";".join(m[6] for m in members if m[6])
                strands = set(m[2] for m in members)
                strand = members[0][2] if len(strands) == 1 else "mixed"
                header = f"{ids or names}|{seqid}:{start}-{end}|{names}|{biotypes}|strand={strand}"
                write_fasta_record(fasta_out, header, seq_str)
                tsv.writerow([seqid, seg_type, start, end, length, strand,
                              names, ids, biotypes, locus_tags, len(members)])
                total_feat_blocks += 1

        if i % 200 == 0 or i == n_seqs:
            print(f"  processed {i}/{n_seqs} sequences...")

    fasta_out.close()
    table_out.close()

    if args.gzip_fasta:
        with open(out_fasta_path, "rb") as fin, gzip.open(out_fasta_path + ".gz", "wb") as fout:
            shutil.copyfileobj(fin, fout)
        os.remove(out_fasta_path)
        out_fasta_path += ".gz"

    print(f"\nDone.")
    print(f"  {args.feature_type} blocks    : {total_feat_blocks}")
    print(f"  {non_feature_label} blocks : {total_other_blocks}")
    print(f"  FASTA : {out_fasta_path}")
    print(f"  Table : {out_table_path}")


if __name__ == "__main__":
    main()