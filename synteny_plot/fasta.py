# -*- coding: utf-8 -*-
"""Sequence record parsing and chromosome/contig filtering."""

from pathlib import Path

from .utils import open_maybe_gzip

FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas"}


def _path_suffixes(path):
    name = Path(path).name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return Path(name).suffixes


def is_fasta_input(path):
    suffixes = _path_suffixes(path)
    return bool(suffixes and suffixes[-1] in FASTA_SUFFIXES)


def _check_records(records, path):
    if not records:
        raise RuntimeError(f"No sequence records found in {path}")

    ids = [r["seq_id"] for r in records]
    dup = sorted(set(x for x in ids if ids.count(x) > 1))
    if dup:
        raise ValueError(f"Duplicated sequence IDs found in {path}: {dup[:20]}")

    return records


def _read_fasta_records(path):
    path = Path(path)

    records = []
    name = None
    length = 0
    order = 0

    with open_maybe_gzip(path, "rt") as f:
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

    return _check_records(records, path)


def _read_fai_records(path):
    path = Path(path)
    records = []
    with open_maybe_gzip(path, "rt") as f:
        for order, line in enumerate(f):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Invalid FAI line in {path}: {line.rstrip()}")
            records.append({"seq_id": parts[0], "length": int(parts[1]), "order": len(records)})
    return _check_records(records, path)


def _read_bed_records(path):
    path = Path(path)
    spans = {}
    order = {}
    with open_maybe_gzip(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Invalid BED line in {path}: {line.rstrip()}")
            chrom = parts[0]
            start = int(float(parts[1]))
            end = int(float(parts[2]))
            if end < start:
                start, end = end, start
            if chrom not in spans:
                spans[chrom] = [start, end]
                order[chrom] = len(order)
            else:
                spans[chrom][0] = min(spans[chrom][0], start)
                spans[chrom][1] = max(spans[chrom][1], end)
    records = [
        {"seq_id": chrom, "length": spans[chrom][1] - spans[chrom][0], "order": order[chrom]}
        for chrom in sorted(order, key=order.get)
    ]
    return _check_records(records, path)


def read_fasta_records(genome_path):
    """Read chromosome records from FASTA, samtools .fai, or BED-like length files."""
    genome_path = Path(genome_path)
    if not genome_path.exists():
        raise FileNotFoundError(f"Genome/record file not found: {genome_path}")

    suffixes = _path_suffixes(genome_path)
    if ".fai" in suffixes:
        return _read_fai_records(genome_path)
    if ".bed" in suffixes:
        return _read_bed_records(genome_path)
    if any(s in FASTA_SUFFIXES for s in suffixes):
        return _read_fasta_records(genome_path)

    with open_maybe_gzip(genome_path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith(">"):
                return _read_fasta_records(genome_path)
            parts = line.split()
            if len(parts) >= 3:
                return _read_bed_records(genome_path)
            if len(parts) >= 2:
                return _read_fai_records(genome_path)
            break

    raise ValueError(f"Cannot detect input format for {genome_path}; use FASTA, .fai, or BED")


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
