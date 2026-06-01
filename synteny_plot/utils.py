# -*- coding: utf-8 -*-
"""Small reusable helpers."""

import gzip
import re
from pathlib import Path

def open_maybe_gzip(path, mode="rt"):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def safe_name(s):
    s = str(s)
    for b in ["/", "\\", " ", "\t", ":", ";", ",", "|", "(", ")", "[", "]", "{", "}"]:
        s = s.replace(b, "_")
    return s


def natural_key(s):
    import re
    parts = re.split(r"(\d+)", str(s))
    return [int(p) if p.isdigit() else p for p in parts]


def split_csv(s):
    if not s:
        return []
    return [x.strip() for x in str(s).split(",") if x.strip()]


def label_from_path(path):
    p = Path(path)
    name = p.name
    for suffix in [".gz", ".fa", ".fasta", ".fna", ".fas"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


def parse_int_list(text):
    if not text:
        return []
    out = []
    for x in split_csv(text):
        out.append(int(float(x)))
    return sorted(set(out))


def parse_float_list(text):
    if not text:
        return []
    out = []
    for x in split_csv(text):
        out.append(float(x))
    return sorted(set(out))
