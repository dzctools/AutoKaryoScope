# -*- coding: utf-8 -*-
"""HTML/CSS/JavaScript template for the interactive SVG plot."""

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8"/>
<title>__TITLE__</title>
<style>
  body {
    margin: 0;
    font-family: Arial, Helvetica, sans-serif;
    background: var(--page-background, #ffffff);
    color: var(--text-color, #111827);
  }
  .toolbar {
    position: sticky;
    top: 0;
    z-index: 10;
    display: flex;
    align-items: center;
    gap: 8px;
    flex-wrap: wrap;
    padding: 10px 14px;
    background: var(--toolbar-background, rgba(255,255,255,0.96));
    border-bottom: 1px solid var(--toolbar-border, #e5e7eb);
    box-shadow: 0 2px 8px rgba(0,0,0,0.035);
  }
  .toolbar label {
    font-size: 13px;
    color: var(--toolbar-label, #374151);
  }
  .toolbar input, .toolbar select {
    border: 1px solid #d1d5db;
    border-radius: 8px;
    padding: 6px 8px;
    font-size: 13px;
    background: white;
  }
  .toolbar button {
    border: 1px solid #cbd5e1;
    border-radius: 8px;
    padding: 6px 12px;
    background: var(--toolbar-button-background, #f8fafc);
    cursor: pointer;
    font-size: 13px;
  }
  .toolbar button:hover {
    background: var(--toolbar-button-hover, #eef2ff);
  }
  .toolbar .check {
    display: inline-flex;
    align-items: center;
    gap: 4px;
  }
  .toolbar .check input {
    padding: 0;
  }
  .svControls {
    display: inline-flex;
    align-items: center;
    gap: 7px;
    flex-wrap: wrap;
  }
  .svSwatch {
    width: 10px;
    height: 10px;
    border-radius: 50%;
    display: inline-block;
    border: 1px solid rgba(17, 24, 39, 0.24);
  }
  .toolbar .sep {
    width: 1px;
    height: 40px;
    background: var(--toolbar-border, #e5e7eb);
    margin: 0 4px;
  }
  #canvasWrap {
    width: 100vw;
    height: calc(100vh - 56px);
    overflow: auto;
    background: var(--canvas-background, white);
  }
  #plotSvg {
    display: block;
    background: var(--canvas-background, white);
  }
  .chr {
    fill: var(--chromosome-color, #1a2847);
    rx: 2;
    ry: 2;
  }
  .chr.odd {
    fill: var(--chromosome-odd-color, #152732);
  }
  .link {
    stroke: none;
    pointer-events: all;
  }
  .link.forward {
    fill: var(--link-forward-color, rgba(180, 180, 180, 0.50));
  }
  .link.reverse {
    fill: var(--link-reverse-color, rgba(214, 120, 95, 0.62));
  }
  .link.thin {
    fill: none;
    pointer-events: stroke;
  }
  .link.forward.thin {
    stroke: var(--link-forward-thin-color, rgba(180, 180, 180, 0.50));
  }
  .link.reverse.thin {
    stroke: var(--link-reverse-thin-color, rgba(214, 120, 95, 0.54));
  }
  .link.svColor {
    stroke: var(--sv-color);
    fill: var(--sv-color);
  }
  .link.svColor.thin {
    fill: none;
  }
  .chrLabel {
    font-size: 12px;
    fill: var(--text-color, #111827);
    text-anchor: middle;
  }
  .scaleTick {
    stroke: var(--scale-tick-color, #6b7280);
    stroke-width: 0.8;
  }
  .scaleLabel {
    font-size: 9px;
    fill: var(--scale-label-color, #6b7280);
    text-anchor: middle;
  }
  .rowLabel {
    font-size: 16px;
    font-weight: 700;
    fill: var(--row-label-color, #111827);
  }
  .titleText {
    font-size: 23px;
    font-weight: 700;
    fill: var(--title-color, #111827);
    text-anchor: middle;
  }
  .subTitle {
    font-size: 12px;
    fill: var(--subtitle-color, #4b5563);
    text-anchor: middle;
  }
  #tooltip {
    position: fixed;
    z-index: 9999;
    pointer-events: none;
    display: none;
    padding: 6px 8px;
    border-radius: 6px;
    background: var(--tooltip-background, rgba(17, 24, 39, 0.92));
    color: var(--tooltip-text, white);
    font-size: 12px;
    line-height: 1.35;
    box-shadow: 0 2px 8px rgba(0,0,0,0.22);
    max-width: 520px;
  }

  #colorPicker {
    position: fixed;
    z-index: 9998;
    display: none;
    padding: 12px;
    border-radius: 10px;
    background: var(--color-picker-background, white);
    box-shadow: 0 4px 20px rgba(0,0,0,0.25);
    border: 1px solid var(--color-picker-border, #e5e7eb);
    width: 340px;
  }
  #colorPicker .cp-title {
    font-size: 12px;
    font-weight: 600;
    margin-bottom: 8px;
    color: #374151;
  }
  #colorPicker .cp-presets {
    display: grid;
    grid-template-columns: repeat(6, 1fr);
    gap: 4px;
    margin-bottom: 10px;
  }
  #colorPicker .cp-swatch {
    width: 40px;
    height: 40px;
    border-radius: 4px;
    cursor: pointer;
    border: 2px solid transparent;
    transition: border-color 0.1s;
  }
  #colorPicker .cp-swatch:hover {
    border-color: #3b82f6;
  }
  #colorPicker .cp-swatch.active {
    border-color: var(--color-picker-active-border, #1d4ed8);
  }
  #colorPicker .cp-row {
    display: flex;
    align-items: center;
    gap: 6px;
    margin-bottom: 6px;
  }
  #colorPicker .cp-row label {
    font-size: 11px;
    color: #6b7280;
    min-width: 36px;
  }
  #colorPicker .cp-row input[type="text"] {
    flex: 1;
    border: 1px solid #d1d5db;
    border-radius: 6px;
    padding: 4px 8px;
    font-size: 12px;
    font-family: monospace;
  }
  #colorPicker .cp-row input[type="color"] {
    width: 30px;
    height: 26px;
    border: none;
    padding: 0;
    cursor: pointer;
  }
  #colorPicker .cp-buttons {
    display: flex;
    gap: 6px;
    margin-top: 8px;
  }
  #colorPicker .cp-buttons button {
    flex: 1;
    border: 1px solid #d1d5db;
    border-radius: 6px;
    padding: 5px 0;
    font-size: 12px;
    cursor: pointer;
    background: #f9fafb;
  }
  #colorPicker .cp-buttons button:hover {
    background: #eef2ff;
  }
  #colorPicker .cp-buttons button.cp-apply {
    background: #3b82f6;
    color: white;
    border-color: #3b82f6;
  }
  #colorPicker .cp-buttons button.cp-apply:hover {
    background: #2563eb;
  }


  #resizeHandle {
    position: fixed;
    right: 0;
    top: 56px;
    bottom: 0;
    width: 14px;
    cursor: ew-resize;
    background: rgba(59, 130, 246, 0.08);
    border-left: 2px solid rgba(59, 130, 246, 0.15);
    z-index: 20;
    transition: background 0.15s, border-color 0.15s;
  }
  #resizeHandle:hover, #resizeHandle.active {
    background: rgba(59, 130, 246, 0.25);
    border-left-color: rgba(59, 130, 246, 0.6);
  }
  #resizeHandle::after {
    content: "⋮⋮";
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%) rotate(90deg);
    font-size: 12px;
    letter-spacing: 2px;
    color: rgba(59, 130, 246, 0.25);
    pointer-events: none;
    transition: color 0.15s;
  }
  #resizeHandle:hover::after, #resizeHandle.active::after {
    color: rgba(59, 130, 246, 0.9);
  }
  @page {
    size: landscape;
    margin: 8mm;
  }
  @media print {
    body {
      background: white;
    }
    .toolbar, #tooltip, .vDragLine {
      display: none !important;
    }
    #canvasWrap {
      overflow: visible !important;
      width: auto !important;
      height: auto !important;
    }
    #plotSvg {
      width: 100% !important;
      height: auto !important;
    }
  }
</style>
</head>
<body>
<div class="toolbar">
  <label>Turn genome</label>
  <select id="turnRow">__ROW_OPTIONS__</select>
  <label>Chromosome</label>
  <input id="turnChr" list="turnChrList" placeholder="HiC_scaffold_1 / chr1" size="34"/>
  <datalist id="turnChrList">__DATALIST__</datalist>
  <button onclick="turnChromosome()">Turn</button>
  <button onclick="resetTurn()">Reset turn</button>

  <span class="sep"></span>

  <label>View genome</label>
  <select id="focusRow">
    <option value="auto">auto</option>
    __ROW_OPTIONS__
  </select>
  <label>Focus chromosome(s)</label>
  <input id="focusChr" list="focusChrList" placeholder="Enter one or more chromosomes, separated by commas" size="40"/>
  <datalist id="focusChrList">__DATALIST__</datalist>
  <button onclick="focusChromosome()">View chr(s)</button>
  <button onclick="resetFocus()">Reset view</button>

  <span class="sep"></span>

  <label class="check"><input id="svColorMode" type="checkbox"/> SV color</label>
  <span id="svTypeControls" class="svControls"></span>

  <span class="sep"></span>

  <button onclick="exportSvg()">Export SVG</button>
  <button onclick="exportPng()">Export PNG</button>
  <button onclick="exportPdf()">Export PDF</button>
  <span class="sep"></span>
  <label>Width:</label>
  <span id="canvasSizeInfo" style="font-size:13px;color:#374151;min-width:60px;"></span>
  <button onclick="resetWidth()">Reset width</button>
</div>

<div id="canvasWrap">
  <svg id="plotSvg" xmlns="http://www.w3.org/2000/svg"></svg>
</div>
<div id="resizeHandle"></div>
<div id="tooltip"></div>
<div id="colorPicker"></div>

<script>
const GENOMES = __GENOMES__;
const BLOCKS = __BLOCKS__;
const SV_EVENTS = __SV_EVENTS__;
const INIT_TURNS = __INIT_TURNS__;
const CONFIG = __CONFIG__;
const COLOR_CONFIG = CONFIG.colors || {};

function cfgColor(key, fallback) {
  return COLOR_CONFIG[key] || fallback;
}

function cfgSvColors() {
  const base = {
    SYN: "#b0b0b0",
    INV: "#dc2626",
    FRAG_INV: "#f97316",
    FUSION: "#dc2626",
    TRANS: "#dc2626",
    INS: "#16a34a",
    DEL: "#f59e0b",
    INS_DEL: "#0d9488",
    OTHER: "#64748b"
  };
  return Object.assign(base, COLOR_CONFIG.sv_colors || COLOR_CONFIG.svColors || {});
}

function applyConfiguredCssColors() {
  const root = document.documentElement;
  const cssMap = {
    "page_background": "--page-background",
    "text_color": "--text-color",
    "toolbar_background": "--toolbar-background",
    "toolbar_label": "--toolbar-label",
    "toolbar_border": "--toolbar-border",
    "toolbar_button_background": "--toolbar-button-background",
    "toolbar_button_hover": "--toolbar-button-hover",
    "canvas_background": "--canvas-background",
    "chromosome": "--chromosome-color",
    "chromosome_odd": "--chromosome-odd-color",
    "link_forward": "--link-forward-color",
    "link_reverse": "--link-reverse-color",
    "link_forward_thin": "--link-forward-thin-color",
    "link_reverse_thin": "--link-reverse-thin-color",
    "scale_tick": "--scale-tick-color",
    "scale_label": "--scale-label-color",
    "row_label": "--row-label-color",
    "title": "--title-color",
    "subtitle": "--subtitle-color",
    "tooltip_background": "--tooltip-background",
    "tooltip_text": "--tooltip-text",
    "color_picker_background": "--color-picker-background",
    "color_picker_border": "--color-picker-border",
    "color_picker_active_border": "--color-picker-active-border"
  };
  for (const [key, cssVar] of Object.entries(cssMap)) {
    if (COLOR_CONFIG[key]) root.style.setProperty(cssVar, COLOR_CONFIG[key]);
  }
}

applyConfiguredCssColors();

let turned = INIT_TURNS.map(x => new Set(x));
let focus = [];  // [{row: number, chrs: ["chr1", ...]}, ...]
let svColorMode = false;

function updateChromDatalist(selectId, datalistId) {
  const rowIdx = parseInt(document.getElementById(selectId).value);
  const datalist = document.getElementById(datalistId);
  if (isNaN(rowIdx) || !GENOMES[rowIdx]) {
    // show all
    const allChr = sortedUniqChr();
    datalist.innerHTML = allChr.map(c => '<option value="' + c + '">').join("\n");
    return;
  }
  const chrs = GENOMES[rowIdx].records.map(r => r.seq_id);
  const sorted = [...new Set(chrs)].sort(naturalSort);
  datalist.innerHTML = sorted.map(c => '<option value="' + c + '">').join("\n");
}

function sortedUniqChr() {
  const all = {};
  for (const g of GENOMES) for (const r of g.records) all[r.seq_id] = 1;
  return Object.keys(all).sort(naturalSort);
}

function naturalSort(a, b) { return a.localeCompare(b, undefined, {numeric: true}); }

document.getElementById("turnRow").addEventListener("change", function() { updateChromDatalist("turnRow", "turnChrList"); });
document.getElementById("focusRow").addEventListener("change", function() { updateChromDatalist("focusRow", "focusChrList"); });
// Init on load
updateChromDatalist("turnRow", "turnChrList");
updateChromDatalist("focusRow", "focusChrList");

let visibleSvTypes = new Set();

const svg = document.getElementById("plotSvg");
const tooltip = document.getElementById("tooltip");

const SV_COLORS = cfgSvColors();

function svTypeOfBlock(b) {
  return String(b.svType || b.LABEL || "SYN").toUpperCase();
}

function svColor(type) {
  return SV_COLORS[type] || SV_COLORS.OTHER;
}

function collectSvTypes() {
  const preferred = ["SYN", "INV", "FRAG_INV", "FUSION", "TRANS", "INS", "DEL", "INS_DEL"];
  const found = new Set(BLOCKS.map(svTypeOfBlock));
  const ordered = preferred.filter(x => found.has(x));
  Array.from(found).sort().forEach(x => {
    if (!ordered.includes(x)) ordered.push(x);
  });
  return ordered.length ? ordered : ["SYN"];
}

function fmtBp(x) {
  x = Number(x);
  if (x >= 1e9) return (x / 1e9).toFixed(2) + " Gb";
  if (x >= 1e6) return (x / 1e6).toFixed(2) + " Mb";
  if (x >= 1e3) return (x / 1e3).toFixed(2) + " kb";
  return Math.round(x) + " bp";
}

function fmtInt(x) {
  return Math.round(Number(x)).toLocaleString("en-US");
}

function invertStrand(s) {
  return s === "+" ? "-" : "+";
}

function clearSvg() {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
}

function makeSvgEl(tag, attrs) {
  const el = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (k === "text") el.textContent = v;
    else el.setAttribute(k, v);
  }
  return el;
}

function showTip(html, event) {
  tooltip.innerHTML = html;
  tooltip.style.display = "block";
  tooltip.style.left = (event.clientX + 12) + "px";
  tooltip.style.top = (event.clientY + 12) + "px";
}

function hideTip() {
  tooltip.style.display = "none";
}

function buildLayout(records, turnedSet, x0, width) {
  if (!records || records.length === 0) return {};
  const totalLen = records.reduce((s, r) => s + r.length, 0);
  const gap = totalLen * CONFIG.gapRatio;
  const total = totalLen + gap * Math.max(0, records.length - 1);

  let cursor = 0;
  const layout = {};

  for (const r of records) {
    const a = x0 + (cursor / total) * width;
    const b = x0 + ((cursor + r.length) / total) * width;

    layout[r.seq_id] = {
      x0: a,
      x1: b,
      len: r.length,
      reversed: turnedSet.has(r.seq_id)
    };

    cursor += r.length + gap;
  }

  return layout;
}

function mapX(layout, chr, bp) {
  const item = layout[chr];
  if (!item) return null;

  bp = Math.max(0, Math.min(bp, item.len));
  let frac = item.len > 0 ? bp / item.len : 0;
  if (item.reversed) frac = 1 - frac;

  return item.x0 + frac * (item.x1 - item.x0);
}

function blockIntervalX(layout, chr, start, end, alnLen) {
  const item = layout[chr];
  if (!item) return null;

  const x1 = mapX(layout, chr, start);
  const x2 = mapX(layout, chr, end);
  if (x1 === null || x2 === null) return null;

  const chrWidth = Math.max(1, Math.abs(item.x1 - item.x0));
  const blockFrac = item.len > 0 ? Math.min(1, Math.max(0, Number(alnLen || 0) / item.len)) : 0;
  const blockWidth = Math.max(1.2, chrWidth * blockFrac);
  const center = (x1 + x2) / 2;
  const left = Math.max(Math.min(item.x0, item.x1), center - blockWidth / 2);
  const right = Math.min(Math.max(item.x0, item.x1), center + blockWidth / 2);

  return {
    left: Math.min(left, right),
    right: Math.max(left, right),
    frac: blockFrac
  };
}

function linkPath(x1, y1, x2, y2, reverse) {
  const midY = (y1 + y2) / 2;
  const bend = CONFIG.curvature * (reverse ? -1 : 1);
  const c1x = x1;
  const c1y = midY + bend;
  const c2x = x2;
  const c2y = midY - bend;
  return `M ${x1.toFixed(2)} ${y1.toFixed(2)} C ${c1x.toFixed(2)} ${c1y.toFixed(2)}, ${c2x.toFixed(2)} ${c2y.toFixed(2)}, ${x2.toFixed(2)} ${y2.toFixed(2)}`;
}

function ribbonPath(u1, u2, y1, l1, l2, y2, reverse) {
  const midY = (y1 + y2) / 2;
  const bend = CONFIG.curvature * (reverse ? -1 : 1);
  const topC1y = midY + bend;
  const topC2y = midY - bend;
  const botC1y = midY - bend;
  const botC2y = midY + bend;
  return (
    `M ${u1.toFixed(2)} ${y1.toFixed(2)} ` +
    `C ${u1.toFixed(2)} ${topC1y.toFixed(2)}, ${l1.toFixed(2)} ${topC2y.toFixed(2)}, ${l1.toFixed(2)} ${y2.toFixed(2)} ` +
    `L ${l2.toFixed(2)} ${y2.toFixed(2)} ` +
    `C ${l2.toFixed(2)} ${botC1y.toFixed(2)}, ${u2.toFixed(2)} ${botC2y.toFixed(2)}, ${u2.toFixed(2)} ${y1.toFixed(2)} Z`
  );
}

function strokeWidth(alnLen) {
  const v = 0.55 + Math.log10(Math.max(alnLen, 1000) / 1000) * 0.42;
  return Math.max(CONFIG.minLineWidth, Math.min(CONFIG.maxLineWidth, v));
}

function drawText(x, y, text, className, anchor="middle") {
  const t = makeSvgEl("text", {x, y, class: className, "text-anchor": anchor});
  const parts = String(text).split("\n");
  parts.forEach((p, i) => {
    const span = makeSvgEl("tspan", {
      x,
      dy: i === 0 ? 0 : 15,
      text: p
    });
    t.appendChild(span);
  });
  svg.appendChild(t);
}


/* ====== Chromosome Color Picker ====== */
const chrColors = {};  // key: "rowIndex:seq_id" -> color string
let cpTarget = null;   // {rowIndex, seq_id, rect}

function chrColorKey(rowIndex, seqId) { return rowIndex + ":" + seqId; }

function showColorPicker(event, rowIndex, seqId) {
  event.stopPropagation();
  const cp = document.getElementById("colorPicker");
  cpTarget = {rowIndex, seqId};
  const key = chrColorKey(rowIndex, seqId);
  const currentColor = chrColors[key] || cfgColor("default_chromosome_color", cfgColor("chromosome", "#1a2847"));

  const presets = COLOR_CONFIG.chromosome_color_presets || COLOR_CONFIG.chromosomeColorPresets || [
    "#1a2847","#152732","#dc2626","#ea580c","#d97706","#ca8a04",
    "#65a30d","#16a34a","#0d9488","#0891b2","#2563eb","#7c3aed",
    "#c026d3","#e11d48","#78716c","#334155","#f8fafc","#1e293b",
  ];

  let html = `<div class="cp-title">${seqId}</div>`;
  html += '<div class="cp-presets">';
  for (const c of presets) {
    const active = c === currentColor ? ' active' : '';
    html += `<div class="cp-swatch${active}" style="background:${c}" data-color="${c}" onclick="selectPresetColor(this,'${c}')"></div>`;
  }
  html += '</div>';
  html += '<div class="cp-row">';
  html += `<input type="color" id="cpNative" value="${currentColor}" onchange="syncColorFromNative()"/>`;
  html += `<input type="text" id="cpInput" value="${currentColor}" placeholder="#dc2626 or 220,38,38" oninput="syncColorFromInput()"/>`;
  html += '</div>';
  html += '<div class="cp-buttons">';
  html += '<button onclick="resetChrColor()">Reset</button>';
  html += '<button onclick="hideColorPicker()">Cancel</button>';
  html += '<button class="cp-apply" onclick="applyChrColor()">Apply</button>';
  html += '</div>';

  cp.innerHTML = html;
  cp.style.display = "block";

  // Position near click
  let left = event.clientX + 12;
  let top = event.clientY + 12;
  if (left + 220 > window.innerWidth) left = event.clientX - 220;
  if (top + 280 > window.innerHeight) top = event.clientY - 280;
  cp.style.left = left + "px";
  cp.style.top = top + "px";
}

function hideColorPicker() {
  document.getElementById("colorPicker").style.display = "none";
  cpTarget = null;
}

function selectPresetColor(el, color) {
  document.querySelectorAll("#colorPicker .cp-swatch").forEach(s => s.classList.remove("active"));
  el.classList.add("active");
  document.getElementById("cpNative").value = color;
  document.getElementById("cpInput").value = color;
}

function syncColorFromNative() {
  document.getElementById("cpInput").value = document.getElementById("cpNative").value;
}

function syncColorFromInput() {
  const val = document.getElementById("cpInput").value.trim();
  // Try HEX
  if (/^#[0-9a-fA-F]{6}$/.test(val)) {
    document.getElementById("cpNative").value = val;
    return;
  }
  // Try RGB: "220,38,38" or "220, 38, 38"
  const m = val.match(/^(\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3})$/);
  if (m) {
    const hex = "#" + [m[1],m[2],m[3]].map(x => parseInt(x).toString(16).padStart(2,"0")).join("");
    document.getElementById("cpNative").value = hex;
  }
}

function hexToRgbStr(hex) {
  const m = hex.match(/^#([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})$/);
  if (!m) return "";
  return `${parseInt(m[1],16)},${parseInt(m[2],16)},${parseInt(m[3],16)}`;
}

function applyChrColor() {
  if (!cpTarget) return;
  const color = document.getElementById("cpInput").value;
  if (/^#[0-9a-fA-F]{6}$/.test(color)) {
    chrColors[chrColorKey(cpTarget.rowIndex, cpTarget.seqId)] = color;
    draw();
  }
  hideColorPicker();
}

function resetChrColor() {
  if (!cpTarget) return;
  delete chrColors[chrColorKey(cpTarget.rowIndex, cpTarget.seqId)];
  draw();
  hideColorPicker();
}

// Close picker on outside click
document.addEventListener("click", function(e) {
  const cp = document.getElementById("colorPicker");
  if (cp.style.display === "block" && !cp.contains(e.target)) {
    hideColorPicker();
  }
});

function drawChromosomes(records, layout, y, rowIndex) {
  const barH = CONFIG.barHeight;
  for (const r of records) {
    const item = layout[r.seq_id];

    const key = chrColorKey(rowIndex, r.seq_id);
    const customColor = chrColors[key] || null;
    const rectAttrs = {
      x: item.x0,
      y: y - barH / 2,
      width: Math.max(1, item.x1 - item.x0),
      height: barH,
      class: "chr",
      style: customColor ? ("fill:" + customColor) : ""
    };
    const rect = makeSvgEl("rect", rectAttrs);

    rect.addEventListener("click", (event) => {
      showColorPicker(event, rowIndex, r.seq_id);
    });
    rect.addEventListener("mousemove", (event) => {
      const pt = svg.createSVGPoint();
      pt.x = event.clientX;
      pt.y = event.clientY;
      const loc = pt.matrixTransform(svg.getScreenCTM().inverse());
      let frac = (loc.x - item.x0) / Math.max(1, item.x1 - item.x0);
      frac = Math.max(0, Math.min(1, frac));
      const bp = item.reversed ? Math.round((1 - frac) * r.length) : Math.round(frac * r.length);
      showTip(`genome: ${GENOMES[rowIndex].label}<br>chr: ${r.seq_id}<br>pos: ${fmtInt(bp)} bp`, event);
    });
    rect.addEventListener("mouseleave", hideTip);

    svg.appendChild(rect);

    const flip = item.reversed ? " ↢" : "";
    const labelY = y - 38;
    // drawText((item.x0 + item.x1) / 2, labelY, `${r.seq_id}${flip}\n${fmtBp(r.length)}`, "chrLabel");

    // Draw scale bar below chromosome
    drawScaleBar(item, r.length, y + barH / 2);
  }
}

function drawScaleBar(item, chrLen, baseY) {
  // Calculate a nice tick interval
  const chrWidth = Math.abs(item.x1 - item.x0);
  if (chrWidth < 30) return; // too small to draw scale

  const bpPerPixel = chrLen / chrWidth;
  // Target: ~3-6 ticks, aim for ~80-150px between ticks
  const targetPx = 100;
  const targetBp = bpPerPixel * targetPx;

  // Round to nice number
  const niceSteps = [1000, 5000, 10000, 20000, 50000, 100000, 200000, 500000,
                     1000000, 2000000, 5000000, 10000000, 20000000, 50000000,
                     100000000, 200000000, 500000000];
  let step = niceSteps[0];
  for (const s of niceSteps) {
    if (s >= targetBp) { step = s; break; }
    step = s;
  }

  const stepPx = (step / chrLen) * chrWidth;
  if (stepPx < 20) return; // ticks too dense

  const numTicks = Math.floor(chrLen / step);
  const tickH = 6;
  const labelOffsetY = tickH + 10;

  for (let i = 0; i <= numTicks; i++) {
    const bp = i * step;
    const frac = bp / chrLen;
    const x = item.x0 + frac * (item.x1 - item.x0);

    // Draw tick line
    const tick = makeSvgEl("line", {
      x1: x, y1: baseY + 2,
      x2: x, y2: baseY + 2 + tickH,
      class: "scaleTick"
    });
    svg.appendChild(tick);

    // Draw label (show for first, last, and some middle ticks to avoid overlap)
    const showLabel = (i === 0 || i === numTicks || stepPx >= 50);
    if (showLabel) {
      const label = makeSvgEl("text", {
        x: x, y: baseY + labelOffsetY,
        class: "scaleLabel",
        "text-anchor": "middle"
      });
      label.textContent = fmtBp(bp);
      svg.appendChild(label);
    }
  }

  // Draw baseline
  const baseline = makeSvgEl("line", {
    x1: item.x0, y1: baseY + 2,
    x2: item.x1, y2: baseY + 2,
    class: "scaleTick"
  });
  svg.appendChild(baseline);
}

function drawLinks(layouts, blocks, rowYByOriginal) {
  const visibleBlocks = blocks.filter(b => visibleSvTypes.has(svTypeOfBlock(b)));
  const wideThreshold = Number(CONFIG.wideRibbonThreshold || 300);
  const pairCounts = new Map();
  for (const b of visibleBlocks) {
    pairCounts.set(b.pair, (pairCounts.get(b.pair) || 0) + 1);
  }
  const sorted = visibleBlocks.slice().sort((a, b) => {
    const aw = (pairCounts.get(a.pair) || 0) < wideThreshold;
    const bw = (pairCounts.get(b.pair) || 0) < wideThreshold;
    if (aw !== bw) return aw ? -1 : 1;
    return aw ? a.alnLen - b.alnLen : b.alnLen - a.alnLen;
  });

  for (const b of sorted) {
    const useWideRibbon = (pairCounts.get(b.pair) || 0) < wideThreshold;
    const upperRow = b.pair;
    const lowerRow = b.pair + 1;

    const upperLayout = layouts[upperRow];
    const lowerLayout = layouts[lowerRow];

    if (!upperLayout || !lowerLayout) continue;
    if (!(b.upper in upperLayout) || !(b.lower in lowerLayout)) continue;

    const upperRev = upperLayout[b.upper].reversed;
    const lowerRev = lowerLayout[b.lower].reversed;

    let visual = b.strand;
    if ((upperRev && !lowerRev) || (!upperRev && lowerRev)) visual = invertStrand(visual);

    const ux1 = mapX(upperLayout, b.upper, b.upperStart);
    const ux2 = mapX(upperLayout, b.upper, b.upperEnd);
    const lx1 = mapX(lowerLayout, b.lower, b.lowerStart);
    const lx2 = mapX(lowerLayout, b.lower, b.lowerEnd);
    if (ux1 === null || ux2 === null || lx1 === null || lx2 === null) continue;

    if (rowYByOriginal[upperRow] === undefined || rowYByOriginal[lowerRow] === undefined) continue;
    const y1 = rowYByOriginal[upperRow] + CONFIG.barHeight / 2;
    const y2 = rowYByOriginal[lowerRow] - CONFIG.barHeight / 2;
    let pathAttrs;
    let upperFillText = "";
    let lowerFillText = "";
    const svType = svTypeOfBlock(b);
    const svStyle = svColorMode ? {"style": `--sv-color: ${svColor(svType)}`} : {};

    if (useWideRibbon) {
      const upperInterval = blockIntervalX(upperLayout, b.upper, b.upperStart, b.upperEnd, b.alnLen);
      const lowerInterval = blockIntervalX(lowerLayout, b.lower, b.lowerStart, b.lowerEnd, b.alnLen);
      if (upperInterval === null || lowerInterval === null) continue;
      const lowerA = visual === "-" ? lowerInterval.right : lowerInterval.left;
      const lowerB = visual === "-" ? lowerInterval.left : lowerInterval.right;
      pathAttrs = {
        d: ribbonPath(upperInterval.left, upperInterval.right, y1, lowerA, lowerB, y2, visual === "-"),
        class: (visual === "+" ? "link forward" : "link reverse") + (svColorMode ? " svColor" : ""),
        "opacity": CONFIG.linkOpacity,
        ...svStyle
      };
      upperFillText = `upper fill: ${(upperInterval.frac * 100).toFixed(2)}%<br>`;
      lowerFillText = `lower fill: ${(lowerInterval.frac * 100).toFixed(2)}%<br>`;
    } else {
      const upperMid = (ux1 + ux2) / 2;
      const lowerMid = (lx1 + lx2) / 2;
      pathAttrs = {
        d: linkPath(upperMid, y1, lowerMid, y2, visual === "-"),
        class: (visual === "+" ? "link forward thin" : "link reverse thin") + (svColorMode ? " svColor" : ""),
        "stroke-width": strokeWidth(b.alnLen),
        "stroke-linecap": "round",
        "opacity": CONFIG.linkOpacity,
        ...svStyle
      };
    }

    const path = makeSvgEl("path", pathAttrs);

    path.addEventListener("mousemove", (event) => {
      showTip(
        `pair: ${GENOMES[upperRow].label} → ${GENOMES[lowerRow].label}<br>` +
        `upper: ${b.upper}:${fmtInt(b.upperStart)}-${fmtInt(b.upperEnd)}<br>` +
        `lower: ${b.lower}:${fmtInt(b.lowerStart)}-${fmtInt(b.lowerEnd)}<br>` +
        `strand: ${b.strand} → displayed ${visual}<br>` +
        `length: ${fmtInt(b.alnLen)} bp<br>` +
        upperFillText +
        lowerFillText +
        `SV: ${svType}${b.svDescription ? " - " + b.svDescription : ""}<br>` +
        `${b.svConfidence !== undefined && b.svConfidence !== "" ? "SV confidence: " + b.svConfidence + "<br>" : ""}` +
        `identity: ${Number(b.identity).toFixed(4)}<br>` +
        `MAPQ: ${b.mapq}`,
        event
      );
    });
    path.addEventListener("mouseleave", hideTip);

    svg.appendChild(path);
  }
}

function parseChrList(text) {
  return String(text || "")
    .split(/[\s,;，；]+/)
    .map(x => x.trim())
    .filter(Boolean);
}

function uniqueKeepOrder(arr) {
  const seen = new Set();
  const out = [];
  for (const x of arr) {
    if (!seen.has(x)) {
      seen.add(x);
      out.push(x);
    }
  }
  return out;
}

function orderRecordsByNames(records, names) {
  const recMap = new Map(records.map(r => [r.seq_id, r]));
  const out = [];
  for (const n of names) {
    if (recMap.has(n)) out.push(recMap.get(n));
  }
  return out;
}

function hasRecord(row, chr) {
  return GENOMES[row].records.some(r => r.seq_id === chr);
}

const FOCUS_PARTNER_MIN_RATIO = 0.20;
const FOCUS_PARTNER_MAX_CHRS = 3;

function selectSupportedFocusPartners(scoreMap) {
  const ranked = Array.from(scoreMap.entries())
    .sort((a, b) => b[1] - a[1] || String(a[0]).localeCompare(String(b[0]), undefined, {numeric: true}));
  if (ranked.length === 0) return [];

  const bestScore = ranked[0][1] || 0;
  const minScore = bestScore * FOCUS_PARTNER_MIN_RATIO;
  return ranked
    .filter(([, score]) => score >= minScore)
    .slice(0, FOCUS_PARTNER_MAX_CHRS)
    .map(([chr]) => chr);
}

function resolveFocusRow(requested, chrList) {
  if (requested !== "auto") return Number(requested);

  let foundRow = null;
  const ambiguous = [];
  const missing = [];

  for (const chr of chrList) {
    const rows = [];
    for (let i = 0; i < GENOMES.length; i++) {
      if (hasRecord(i, chr)) rows.push(i);
    }
    if (rows.length === 0) missing.push(chr);
    else if (rows.length > 1) ambiguous.push(chr);
    else if (foundRow === null) foundRow = rows[0];
    else if (foundRow !== rows[0]) ambiguous.push(chr);
  }

  if (missing.length > 0) {
    alert("These chromosomes were not found:\n" + missing.join(", "));
    return null;
  }
  if (ambiguous.length > 0) {
    alert("Auto detection failed. These chromosomes exist in multiple genomes, or the input mixes chromosomes from different genomes. Please select View genome manually.\n" + ambiguous.join(", "));
    return null;
  }

  return foundRow;
}

function getCurrentData() {
  if (!focus || focus.length === 0) {
    return {
      genomes: GENOMES.map((g, i) => ({label: g.label, records: g.records, rowIndex: i})),
      blocks: BLOCKS,
      note: ""
    };
  }

  const n = GENOMES.length;
  const active = [];
  for (let i = 0; i < n; i++) active.push(new Set());

  for (const item of focus) {
    const selectedRow = item.row;
    let upperRow = selectedRow;
    let lowerRow = selectedRow + 1;
    if (selectedRow >= n - 1) {
      upperRow = Math.max(0, selectedRow - 1);
      lowerRow = selectedRow;
    }

    const selectedChrs = new Set(item.chrs);
    for (const chr of item.chrs) active[selectedRow].add(chr);

    const partnerScores = new Map();
    for (const b of BLOCKS) {
      if (b.pair !== upperRow) continue;
      if (selectedRow === upperRow && selectedChrs.has(b.upper)) {
        partnerScores.set(b.lower, (partnerScores.get(b.lower) || 0) + Number(b.alnLen || 0));
      }
      if (selectedRow === lowerRow && selectedChrs.has(b.lower)) {
        partnerScores.set(b.upper, (partnerScores.get(b.upper) || 0) + Number(b.alnLen || 0));
      }
    }

    const selectedPartners = selectSupportedFocusPartners(partnerScores);
    const partnerRow = selectedRow === upperRow ? lowerRow : upperRow;
    for (const chr of selectedPartners) active[partnerRow].add(chr);
  }

  const genomes = [];
  for (let i = 0; i < n; i++) {
    if (active[i].size === 0) continue;
    const names = Array.from(active[i]);
    genomes.push({
      label: GENOMES[i].label,
      records: orderRecordsByNames(GENOMES[i].records, names),
      rowIndex: i
    });
  }

  const blocks = BLOCKS.filter(b => {
    const upperRow = b.pair;
    const lowerRow = b.pair + 1;
    return active[upperRow] &&
      active[lowerRow] &&
      active[upperRow].has(b.upper) &&
      active[lowerRow].has(b.lower);
  });

  const focusText = focus
    .map(item => `${GENOMES[item.row].label} [${item.chrs.join(", ")}]`)
    .join("; ");

  return {
    genomes,
    blocks,
    note: `Focused: ${focusText}, ${blocks.length} block(s)`
  };
}

function setupSvControls() {
  const types = collectSvTypes();
  visibleSvTypes = new Set(types);

  const wrap = document.getElementById("svTypeControls");
  wrap.innerHTML = "";

  for (const type of types) {
    const label = document.createElement("label");
    label.className = "check";

    const input = document.createElement("input");
    input.type = "checkbox";
    input.checked = true;
    input.dataset.svType = type;
    input.addEventListener("change", () => {
      if (input.checked) visibleSvTypes.add(type);
      else visibleSvTypes.delete(type);
      draw();
    });

    const swatch = document.createElement("span");
    swatch.className = "svSwatch";
    swatch.style.background = svColor(type);

    const text = document.createElement("span");
    text.textContent = type;

    label.appendChild(input);
    label.appendChild(swatch);
    label.appendChild(text);
    wrap.appendChild(label);
  }

  document.getElementById("svColorMode").addEventListener("change", (event) => {
    svColorMode = event.target.checked;
    draw();
  });
}

function draw() {
  clearSvg();

  svg.setAttribute("width", CONFIG.width);
  svg.setAttribute("height", CONFIG.height);
  svg.setAttribute("viewBox", `0 0 ${CONFIG.width} ${CONFIG.height}`);

  const current = getCurrentData();

  drawText(CONFIG.width / 2, 34, CONFIG.title, "titleText");
  if (current.note) {
    drawText(CONFIG.width / 2, 58, current.note, "subTitle");
  }

  const plotWidth = CONFIG.width - CONFIG.marginLeft - CONFIG.marginRight;
  const layouts = [];
  const rowYByOriginal = [];

  for (let i = 0; i < current.genomes.length; i++) {
    const originalRow = current.genomes[i].rowIndex !== undefined ? current.genomes[i].rowIndex : i;
    const y = CONFIG.rowYs[i];
    drawText(CONFIG.marginLeft - 18, y + 5, current.genomes[i].label, "rowLabel", "end");
    const layout = buildLayout(current.genomes[i].records, turned[originalRow], CONFIG.marginLeft, plotWidth);
    layouts[originalRow] = layout;
    rowYByOriginal[originalRow] = y;
  }

  drawLinks(layouts, current.blocks, rowYByOriginal);

  for (let i = 0; i < current.genomes.length; i++) {
    const originalRow = current.genomes[i].rowIndex !== undefined ? current.genomes[i].rowIndex : i;
    drawChromosomes(current.genomes[i].records, layouts[originalRow], CONFIG.rowYs[i], originalRow);
  }

  const visibleCount = current.blocks.filter(b => visibleSvTypes.has(svTypeOfBlock(b))).length;
  if (visibleCount === 0) {
    drawText(CONFIG.width / 2, CONFIG.height / 2, "No synteny blocks under current filters/view", "subTitle");
  }

  // Draw draggable green lines between rows
  for (let i = 0; i < current.genomes.length - 1; i++) {
    const y1 = CONFIG.rowYs[i];
    const y2 = CONFIG.rowYs[i + 1];
    const midY = (y1 + y2) / 2;

    const g = makeSvgEl("g", {class: "vDragLine"});

    // Hit area (invisible, wider for easy grabbing)
    const hitArea = makeSvgEl("rect", {
      x: 0, y: midY - 12,
      width: CONFIG.width, height: 24,
      fill: "transparent",
      cursor: "ns-resize",
      style: "pointer-events: all;"
    });

    // Visible green line
    const line = makeSvgEl("line", {
      x1: CONFIG.marginLeft, y1: midY,
      x2: CONFIG.width - CONFIG.marginRight, y2: midY,
      stroke: "rgba(16, 185, 129, 0.4)",
      "stroke-width": 2,
      "stroke-dasharray": "8,4",
      style: "pointer-events: none;"
    });

    // Label
    const label = makeSvgEl("text", {
      x: CONFIG.marginLeft - 8, y: midY + 4,
      fill: "rgba(16, 185, 129, 0.6)",
      "font-size": "11px",
      "text-anchor": "end",
      style: "pointer-events: none;"
    });
    label.textContent = "⇕ drag";

    g.appendChild(hitArea);
    g.appendChild(line);
    g.appendChild(label);

    // Drag logic
    let startY, startYs;
    hitArea.addEventListener("mousedown", (e) => {
      e.preventDefault();
      startY = e.clientY;
      startYs = CONFIG.rowYs.slice();
      line.setAttribute("stroke", "rgba(16, 185, 129, 0.9)");
      line.setAttribute("stroke-width", "3");
      const onMove = (ev) => {
        const dy = ev.clientY - startY;
        const gap = startYs[i + 1] - startYs[i];
        const newGap = Math.max(80, gap + dy);
        const actualDy = newGap - gap;
        for (let j = i + 1; j < CONFIG.rowYs.length; j++) {
          CONFIG.rowYs[j] = startYs[j] + actualDy;
        }
        CONFIG.height = Math.max(CONFIG.height, CONFIG.rowYs[CONFIG.rowYs.length - 1] + 100);
        draw();
      };
      const onUp = () => {
        document.removeEventListener("mousemove", onMove);
        document.removeEventListener("mouseup", onUp);
        document.body.style.cursor = "";
        document.body.style.userSelect = "";
      };
      document.body.style.cursor = "ns-resize";
      document.body.style.userSelect = "none";
      document.addEventListener("mousemove", onMove);
      document.addEventListener("mouseup", onUp);
    });

    svg.appendChild(g);
  }
}


/* ====== Drag-to-resize canvas ====== */
(function() {
  const handle = document.getElementById("resizeHandle");
  const ORIGINAL_WIDTH = CONFIG.width;
  const ORIGINAL_ROW_YS = CONFIG.rowYs.slice();
  let dragging = false;
  let startX, startWidth;

  function onDown(e) {
    e.preventDefault();
    dragging = true;
    startX = e.clientX || (e.touches && e.touches[0].clientX);
    startWidth = CONFIG.width;
    handle.classList.add("active");
    document.body.style.cursor = "ew-resize";
    document.body.style.userSelect = "none";
  }

  function onMove(e) {
    if (!dragging) return;
    const clientX = e.clientX || (e.touches && e.touches[0].clientX);
    const dx = clientX - startX;
    const newWidth = Math.max(600, startWidth + dx);
    CONFIG.width = newWidth;
    draw();
    const sizeInfo = document.getElementById("canvasSizeInfo");
    if (sizeInfo) sizeInfo.textContent = newWidth + "px";
  }

  function onUp() {
    if (!dragging) return;
    dragging = false;
    handle.classList.remove("active");
    document.body.style.cursor = "";
    document.body.style.userSelect = "";
  }

  handle.addEventListener("mousedown", onDown);
  document.addEventListener("mousemove", onMove);
  document.addEventListener("mouseup", onUp);
  handle.addEventListener("touchstart", onDown, {passive: false});
  document.addEventListener("touchmove", onMove, {passive: false});
  document.addEventListener("touchend", onUp);
})();

function turnChromosome() {
  const row = Number(document.getElementById("turnRow").value);
  const chr = document.getElementById("turnChr").value.trim();

  if (!chr) {
    alert("Please enter the chromosome name to turn.");
    return;
  }
  if (!hasRecord(row, chr)) {
    alert(`Not found in ${GENOMES[row].label}: ${chr}`);
    return;
  }

  if (turned[row].has(chr)) turned[row].delete(chr);
  else turned[row].add(chr);

  draw();
}

function resetTurn() {
  turned = INIT_TURNS.map(x => new Set(x));
  draw();
}

function focusChromosome() {
  const requested = document.getElementById("focusRow").value;
  const raw = document.getElementById("focusChr").value.trim();

  if (!raw) {
    alert("Please enter one or more chromosome names to view, separated by commas.");
    return;
  }

  const chrList = uniqueKeepOrder(parseChrList(raw));
  const row = resolveFocusRow(requested, chrList);
  if (row === null || row === undefined || Number.isNaN(row)) return;

  const missing = chrList.filter(chr => !hasRecord(row, chr));
  if (missing.length > 0) {
    alert(`Not found in ${GENOMES[row].label}:\n` + missing.join(", "));
    return;
  }

  const existingFocus = focus.find(item => item.row === row);
  const mergedChrList = existingFocus
    ? uniqueKeepOrder(existingFocus.chrs.concat(chrList))
    : chrList;
  if (existingFocus) {
    existingFocus.chrs = mergedChrList;
  } else {
    focus.push({row, chrs: mergedChrList});
  }
  document.getElementById("focusRow").value = String(row);
  document.getElementById("focusChr").value = mergedChrList.join(", ");
  draw();
}

function resetFocus() {
  focus = [];
  document.getElementById("focusChr").value = "";
  document.getElementById("focusRow").value = "auto";
  draw();
}

function exportPdf() {
  hideTip();
  // Browser will open the print dialog. Choose "Save as PDF".
  window.print();
}

function exportFileName(ext) {
  const safeTitle = String(CONFIG.title || "multi_synteny")
    .replace(/[\\/:*?"<>|]+/g, "_")
    .replace(/\s+/g, "_")
    .replace(/^_+|_+$/g, "") || "multi_synteny";
  return `${safeTitle}.${ext}`;
}

function cloneSvgForExport() {
  const clone = svg.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("width", CONFIG.width);
  clone.setAttribute("height", CONFIG.height);
  clone.setAttribute("viewBox", `0 0 ${CONFIG.width} ${CONFIG.height}`);

  inlineComputedSvgStyles(svg, clone);
  removeExportOnlyHelpers(clone);

  const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
  style.textContent = Array.from(document.styleSheets)
    .map(sheet => {
      try {
        return Array.from(sheet.cssRules).map(rule => rule.cssText).join("\n");
      } catch (err) {
        return "";
      }
    })
    .filter(Boolean)
    .join("\n");
  clone.insertBefore(style, clone.firstChild);

  const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  bg.setAttribute("x", 0);
  bg.setAttribute("y", 0);
  bg.setAttribute("width", CONFIG.width);
  bg.setAttribute("height", CONFIG.height);
  bg.setAttribute("fill", "white");
  clone.insertBefore(bg, clone.firstChild);
  return clone;
}

function removeExportOnlyHelpers(cloneRoot) {
  cloneRoot.querySelectorAll(".vDragLine").forEach(el => el.remove());
}

function inlineComputedSvgStyles(sourceRoot, cloneRoot) {
  const attrProps = [
    "stroke-width", "stroke-linecap", "stroke-linejoin",
    "font-family", "font-size", "font-weight", "text-anchor"
  ];

  function copyStyles(sourceNode, cloneNode) {
    if (sourceNode.nodeType === 1 && cloneNode.nodeType === 1) {
      const computed = window.getComputedStyle(sourceNode);

      setSvgColorAttr(cloneNode, "fill", computed.getPropertyValue("fill"), computed.getPropertyValue("fill-opacity"));
      setSvgColorAttr(cloneNode, "stroke", computed.getPropertyValue("stroke"), computed.getPropertyValue("stroke-opacity"));

      const opacity = computed.getPropertyValue("opacity");
      if (opacity && opacity !== "1") cloneNode.setAttribute("opacity", opacity);

      for (const prop of attrProps) {
        const value = computed.getPropertyValue(prop);
        if (value && value !== "auto" && value !== "normal") {
          cloneNode.setAttribute(prop, value);
        }
      }

      cloneNode.removeAttribute("class");
      cloneNode.removeAttribute("style");
    }

    const sourceChildren = sourceNode.childNodes;
    const cloneChildren = cloneNode.childNodes;
    for (let i = 0; i < sourceChildren.length; i++) {
      copyStyles(sourceChildren[i], cloneChildren[i]);
    }
  }

  copyStyles(sourceRoot, cloneRoot);
}

function setSvgColorAttr(el, attr, colorValue, opacityValue) {
  if (!colorValue) return;
  const normalized = normalizeSvgColor(colorValue);
  if (!normalized) return;

  el.setAttribute(attr, normalized.color);
  const combinedOpacity = Number(normalized.opacity) * Number(opacityValue || 1);
  if (combinedOpacity < 1) {
    el.setAttribute(`${attr}-opacity`, String(Math.max(0, Math.min(1, combinedOpacity))));
  }
}

function normalizeSvgColor(value) {
  const v = String(value).trim();
  if (!v || v === "none" || v === "transparent") return {color: "none", opacity: 1};

  const rgba = v.match(/^rgba?\(([^)]+)\)$/i);
  if (!rgba) return {color: v, opacity: 1};

  const parts = rgba[1].split(",").map(x => x.trim());
  if (parts.length < 3) return {color: v, opacity: 1};

  const rgb = parts.slice(0, 3).map(x => {
    if (x.endsWith("%")) return Math.round(parseFloat(x) * 2.55);
    return Math.round(parseFloat(x));
  }).map(x => Math.max(0, Math.min(255, x)));
  const alpha = parts.length >= 4 ? parseFloat(parts[3]) : 1;
  const hex = "#" + rgb.map(x => x.toString(16).padStart(2, "0")).join("");
  return {color: hex, opacity: Number.isFinite(alpha) ? alpha : 1};
}

function serializedSvg() {
  const clone = cloneSvgForExport();
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
}

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function exportSvg() {
  hideTip();
  const blob = new Blob([serializedSvg()], {type: "image/svg+xml;charset=utf-8"});
  downloadBlob(blob, exportFileName("svg"));
}

function exportPng() {
  hideTip();
  const svgText = serializedSvg();
  const svgBlob = new Blob([svgText], {type: "image/svg+xml;charset=utf-8"});
  const url = URL.createObjectURL(svgBlob);
  const img = new Image();
  img.onload = function() {
    const canvas = document.createElement("canvas");
    canvas.width = CONFIG.width;
    canvas.height = CONFIG.height;
    const ctx = canvas.getContext("2d");
    ctx.fillStyle = "white";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0);
    canvas.toBlob(blob => {
      URL.revokeObjectURL(url);
      if (blob) downloadBlob(blob, exportFileName("png"));
    }, "image/png");
  };
  img.onerror = function() {
    URL.revokeObjectURL(url);
    alert("PNG export failed. Please try SVG export.");
  };
  img.src = url;
}

document.getElementById("turnChr").addEventListener("keydown", function(e) {
  if (e.key === "Enter") turnChromosome();
});

document.getElementById("focusChr").addEventListener("keydown", function(e) {
  if (e.key === "Enter") focusChromosome();
});

setupSvControls();
draw();
</script>
</body>
</html>
"""
