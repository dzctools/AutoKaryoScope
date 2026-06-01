# -*- coding: utf-8 -*-
"""User configuration loading for command defaults and HTML colors."""

import json
from copy import deepcopy
from pathlib import Path


DEFAULT_COLORS = {
    "page_background": "#ffffff",
    "text_color": "#111827",
    "toolbar_background": "rgba(255,255,255,0.96)",
    "toolbar_label": "#374151",
    "toolbar_border": "#e5e7eb",
    "toolbar_button_background": "#f8fafc",
    "toolbar_button_hover": "#eef2ff",
    "canvas_background": "#ffffff",
    "chromosome": "#1a2847",
    "chromosome_odd": "#152732",
    "link_forward": "rgba(180, 180, 180, 0.50)",
    "link_reverse": "rgba(214, 120, 95, 0.62)",
    "link_forward_thin": "rgba(180, 180, 180, 0.50)",
    "link_reverse_thin": "rgba(214, 120, 95, 0.54)",
    "scale_tick": "#6b7280",
    "scale_label": "#6b7280",
    "row_label": "#111827",
    "title": "#111827",
    "subtitle": "#4b5563",
    "tooltip_background": "rgba(17, 24, 39, 0.92)",
    "tooltip_text": "#ffffff",
    "color_picker_background": "#ffffff",
    "color_picker_border": "#e5e7eb",
    "color_picker_active_border": "#1d4ed8",
    "default_chromosome_color": "#1a2847",
    "chromosome_color_presets": [
        "#1a2847", "#152732", "#dc2626", "#ea580c", "#d97706", "#ca8a04",
        "#65a30d", "#16a34a", "#0d9488", "#0891b2", "#2563eb", "#7c3aed",
        "#c026d3", "#e11d48", "#78716c", "#334155", "#f8fafc", "#1e293b",
    ],
    "sv_colors": {
        "SYN": "#b0b0b0",
        "INV": "#dc2626",
        "FRAG_INV": "#f97316",
        "FUSION": "#dc2626",
        "TRANS": "#dc2626",
        "INS": "#16a34a",
        "DEL": "#f59e0b",
        "INS_DEL": "#0d9488",
        "OTHER": "#64748b",
    },
}


CONFIG_TEMPLATE = {
    "parameters": {
        "width": 3500,
        "height": 0,
        "margin_left": None,
        "margin_right": 70,
        "top_y": 150,
        "row_gap": 260,
        "bar_height": 24,
        "gap_ratio": 0.012,
        "curvature": 80.0,
        "link_opacity": 0.58,
        "min_line_width": 0.45,
        "max_line_width": 3.0,
        "wide_ribbon_threshold": 300,
        "min_chr_len": 3000000,
        "min_block": 50000,
        "min_mapq": 60,
        "min_identity": 0.90,
        "min_pair_aln": 10000000,
        "max_partners_per_chr": 4,
        "partner_cap_mode": "both",
        "final_noise_filter": True,
        "final_noise_continuity_gap": 5000000,
        "final_noise_borderline_ratio": 1.25,
    },
    "colors": DEFAULT_COLORS,
}


def deep_merge(base, override):
    result = deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_user_config(path):
    if not path:
        return {"parameters": {}, "colors": deepcopy(DEFAULT_COLORS)}

    cfg_path = Path(path).expanduser().resolve()
    data = json.loads(cfg_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Config must be a JSON object: {cfg_path}")

    parameters = data.get("parameters", data.get("params", data.get("defaults", {}))) or {}
    colors = deep_merge(DEFAULT_COLORS, data.get("colors", {}) or {})
    if not isinstance(parameters, dict):
        raise ValueError("Config field 'parameters' must be a JSON object")
    if not isinstance(colors, dict):
        raise ValueError("Config field 'colors' must be a JSON object")
    return {"parameters": parameters, "colors": colors}


def write_config_template(path):
    cfg_path = Path(path).expanduser().resolve()
    cfg_path.write_text(
        json.dumps(CONFIG_TEMPLATE, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return cfg_path
