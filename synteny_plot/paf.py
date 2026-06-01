# -*- coding: utf-8 -*-
"""Minimap2 execution and PAF parsing."""

import shutil
import subprocess
import sys
from pathlib import Path

def run_minimap2(target_fa, query_fa, paf_out,
                 minimap2="minimap2",
                 threads=16,
                 preset="asm5",
                 force=False):
    paf_out = Path(paf_out)

    if paf_out.exists() and paf_out.stat().st_size > 0 and not force:
        print(f"[reuse] PAF exists: {paf_out}", file=sys.stderr)
        return

    if shutil.which(minimap2) is None:
        raise RuntimeError(
            f"Cannot find minimap2: {minimap2}. "
            "Load minimap2 into PATH or provide existing PAFs with --pafs."
        )

    cmd = [
        minimap2,
        "-x", preset,
        "--secondary=no",
        "--cs",
        "-c",
        "--eqx",
        "-t", str(threads),
        str(target_fa),
        str(query_fa),
    ]

    print("[run]", " ".join(cmd), ">", paf_out, file=sys.stderr)
    with paf_out.open("w") as out:
        subprocess.run(cmd, stdout=out, stderr=sys.stderr, check=True)


def parse_paf_for_pair(paf_file, upper_ids, lower_ids,
                       pair_index,
                       upper_label,
                       lower_label,
                       min_block=10000,
                       min_mapq=0,
                       min_identity=0.65):
    """
    PAF:
      normal orientation:
        query  = lower genome of this pair
        target = upper genome of this pair
      reversed orientation is also accepted:
        query  = upper genome
        target = lower genome

    Returned block fields are generalized:
      pair, upper, lower, upperStart, lowerStart, ...
    """
    upper_ids = set(upper_ids)
    lower_ids = set(lower_ids)
    blocks = []

    with Path(paf_file).open() as f:
        for line in f:
            if not line.strip():
                continue

            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue

            qname = p[0]
            qlen = int(p[1])
            qstart = int(p[2])
            qend = int(p[3])
            strand = p[4]
            tname = p[5]
            tlen = int(p[6])
            tstart = int(p[7])
            tend = int(p[8])
            nmatch = int(p[9])
            alnlen = int(p[10])
            mapq = int(p[11])

            is_normal = tname in upper_ids and qname in lower_ids
            is_reversed = qname in upper_ids and tname in lower_ids
            if not is_normal and not is_reversed:
                continue
            if alnlen < min_block:
                continue
            if mapq < min_mapq:
                continue

            identity = nmatch / alnlen if alnlen else 0.0
            if identity < min_identity:
                continue

            if is_normal:
                upper = tname
                upper_start = tstart
                upper_end = tend
                upper_len = tlen
                lower = qname
                lower_start = qstart
                lower_end = qend
                lower_len = qlen
            else:
                upper = qname
                upper_start = qstart
                upper_end = qend
                upper_len = qlen
                lower = tname
                lower_start = tstart
                lower_end = tend
                lower_len = tlen

            blocks.append({
                "pair": pair_index,
                "upperGenome": upper_label,
                "lowerGenome": lower_label,
                "upper": upper,
                "upperStart": upper_start,
                "upperEnd": upper_end,
                "upperLen": upper_len,
                "lower": lower,
                "lowerStart": lower_start,
                "lowerEnd": lower_end,
                "lowerLen": lower_len,
                "strand": strand,
                "alnLen": alnlen,
                "identity": round(identity, 5),
                "mapq": mapq,
            })

    return blocks
