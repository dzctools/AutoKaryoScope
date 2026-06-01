# -*- coding: utf-8 -*-
"""FASTA parsing and chromosome/contig filtering."""

from pathlib import Path

from .utils import open_maybe_gzip

def read_fasta_records(fasta_path):
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    records = []
    name = None
    length = 0
    order = 0

    with open_maybe_gzip(fasta_path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records.append({"seq_id": name, "length": length, "order": order})
                    order += 1
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())

    if name is not None:
        records.append({"seq_id": name, "length": length, "order": order})

    if not records:
        raise RuntimeError(f"No FASTA records found in {fasta_path}")

    ids = [r["seq_id"] for r in records]
    dup = sorted(set(x for x in ids if ids.count(x) > 1))
    if dup:
        raise ValueError(f"Duplicated FASTA IDs found in {fasta_path}: {dup[:20]}")

    return records


def filter_records_by_length(records, min_chr_len):
    return [r for r in records if r["length"] >= min_chr_len]


def auto_min_chr_len_for_records(records, floor=2_000_000, ratio=0.01):
    """
    Automatically choose a chromosome-like contig length threshold.

    Default rule:
      max(floor, longest_contig * ratio)

    This removes small contigs/scaffolds before PAF parsing and plotting.
    """
    if not records:
        return int(floor)
    longest = max(int(r["length"]) for r in records)
    return int(max(floor, longest * ratio))
