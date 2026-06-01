# -*- coding: utf-8 -*-
"""HTML assembly from genomes, blocks, turns, SV events, and display configuration."""

import json

from .html_template import HTML_TEMPLATE
from .utils import natural_key

def build_html(genomes, blocks, init_turns,
               title,
               width, height,
               margin_left, margin_right,
               row_ys,
               bar_height,
               gap_ratio,
               curvature,
               link_opacity,
               min_line_width,
               max_line_width,
               wide_ribbon_threshold=300,
               sv_events=None,
               colors=None):
    all_chr = sorted({r["seq_id"] for g in genomes for r in g["records"]}, key=natural_key)
    datalist = "\n".join(f'<option value="{x}"></option>' for x in all_chr)

    row_options = "\n".join(
        f'<option value="{i}">{i}: {g["label"]}</option>'
        for i, g in enumerate(genomes)
    )

    config = {
        "title": title,
        "width": width,
        "height": height,
        "marginLeft": margin_left,
        "marginRight": margin_right,
        "rowYs": row_ys,
        "barHeight": bar_height,
        "gapRatio": gap_ratio,
        "curvature": curvature,
        "linkOpacity": link_opacity,
        "minLineWidth": min_line_width,
        "maxLineWidth": max_line_width,
        "wideRibbonThreshold": int(wide_ribbon_threshold),
        "colors": colors or {},
    }

    sv_data = [e.to_dict() for e in (sv_events or [])]

    html = HTML_TEMPLATE
    html = html.replace("__TITLE__", title)
    html = html.replace("__ROW_OPTIONS__", row_options)
    html = html.replace("__DATALIST__", datalist)
    html = html.replace("__GENOMES__", json.dumps(genomes, ensure_ascii=False))
    html = html.replace("__BLOCKS__", json.dumps(blocks, ensure_ascii=False))
    html = html.replace("__INIT_TURNS__", json.dumps([sorted(x) for x in init_turns], ensure_ascii=False))
    html = html.replace("__CONFIG__", json.dumps(config, ensure_ascii=False))
    html = html.replace("__SV_EVENTS__", json.dumps(sv_data, ensure_ascii=False))
    return html
