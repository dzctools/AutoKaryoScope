# -*- coding: utf-8 -*-
"""Data structures for structural variant events."""


class SVEvent:
    """Represents a detected structural variant between two genomes."""

    __slots__ = (
        "sv_type", "pair",
        "ref_chr", "ref_pos",
        "qry_chr", "qry_pos",
        "size", "confidence",
        "support_blocks", "description",
    )

    def __init__(
        self,
        sv_type,
        pair,
        ref_chr,
        ref_pos,
        qry_chr,
        qry_pos,
        size,
        confidence=0.0,
        support_blocks=None,
        description="",
    ):
        self.sv_type = sv_type
        self.pair = pair
        self.ref_chr = ref_chr
        self.ref_pos = ref_pos
        self.qry_chr = qry_chr
        self.qry_pos = qry_pos
        self.size = size
        self.confidence = confidence
        self.support_blocks = support_blocks or []
        self.description = description

    def to_dict(self):
        return {
            "sv_type": self.sv_type,
            "pair": self.pair,
            "ref_chr": self.ref_chr,
            "ref_pos": self.ref_pos,
            "qry_chr": self.qry_chr,
            "qry_pos": self.qry_pos,
            "size": self.size,
            "confidence": round(self.confidence, 4),
            "description": self.description,
            "block_count": len(self.support_blocks),
        }
