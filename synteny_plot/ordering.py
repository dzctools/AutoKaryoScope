# -*- coding: utf-8 -*-
"""Chromosome ordering and display-only flip parsing."""

from .filters import count_chr_pair_blocks
from .utils import natural_key, split_csv

def choose_order(records, matched_ids, mode="input"):
    ids = [r["seq_id"] for r in records if r["seq_id"] in matched_ids]
    rec = {r["seq_id"]: r for r in records}

    if mode == "input":
        return ids
    if mode == "name":
        return sorted(ids, key=natural_key)
    if mode == "length":
        return sorted(ids, key=lambda x: rec[x]["length"], reverse=True)

    raise ValueError(f"Unknown order mode: {mode}")


def choose_synteny_orders(records_by_row, matched_by_row, blocks):
    """
    Order chromosomes so dominant matched partners appear as vertically as possible.

    Row 0 keeps FASTA/input order. Each lower row is sorted according to the position of its
    strongest upper-row partner in the already-ordered row above.
    """
    n = len(records_by_row)
    if n == 0:
        return []

    input_index = [
        {r["seq_id"]: idx for idx, r in enumerate(records_by_row[i])}
        for i in range(n)
    ]

    row_orders = []
    row0 = [r["seq_id"] for r in records_by_row[0] if r["seq_id"] in matched_by_row[0]]
    row_orders.append(row0)

    for row in range(1, n):
        ids = [r["seq_id"] for r in records_by_row[row] if r["seq_id"] in matched_by_row[row]]
        prev_pos = {c: i for i, c in enumerate(row_orders[row - 1])}

        pair_blocks = [b for b in blocks if int(b["pair"]) == row - 1]
        stat = count_chr_pair_blocks(pair_blocks)

        best_upper_for_lower = {}
        for x in stat.values():
            lower = x["lower"]
            candidate = (x["upper"], x["aln_sum"], x["block_count"])
            if lower not in best_upper_for_lower:
                best_upper_for_lower[lower] = candidate
            else:
                old = best_upper_for_lower[lower]
                if (candidate[1], candidate[2]) > (old[1], old[2]):
                    best_upper_for_lower[lower] = candidate

        def sort_key(chr_id):
            if chr_id in best_upper_for_lower:
                upper, aln_sum, block_count = best_upper_for_lower[chr_id]
                return (
                    0,
                    prev_pos.get(upper, 10**9),
                    -int(aln_sum),
                    -int(block_count),
                    natural_key(chr_id),
                )
            return (1, input_index[row].get(chr_id, 10**9), 0, 0, natural_key(chr_id))

        row_orders.append(sorted(ids, key=sort_key))

    return row_orders


def parse_turn_spec_multi(turn, labels, n_genomes):
    """
    Supported:
      --turn chr1,chr2
          flips chromosomes in genome row 0

      --turn '0:chr1,chr2;2:chr5'
          flips chr1/chr2 in row 0 and chr5 in row 2

      --turn 'top:chr1;bottom:chr2'
          aliases for row 0 and row 1

      --turn 'T2T:chr1;Pub:chr2'
          uses genome label as row selector
    """
    out = {i: set() for i in range(n_genomes)}
    turn = (turn or "").strip()

    if not turn or turn.lower() in {"none", "na", "null", "-"}:
        return out

    def row_from_key(key):
        key = key.strip()
        if key == "top":
            return 0
        if key == "bottom":
            return 1 if n_genomes > 1 else None
        if key.isdigit():
            i = int(key)
            return i if 0 <= i < n_genomes else None
        for i, lab in enumerate(labels):
            if key == lab:
                return i
        return None

    if ":" not in turn:
        for x in split_csv(turn):
            out[0].add(x)
        return out

    for part in [x.strip() for x in turn.split(";") if x.strip()]:
        if ":" not in part:
            for x in split_csv(part):
                out[0].add(x)
            continue

        key, vals = part.split(":", 1)
        row = row_from_key(key)
        if row is None:
            print(f"[warning] cannot resolve --turn group: {key}", file=sys.stderr)
            continue
        for x in split_csv(vals):
            out[row].add(x)

    return out
