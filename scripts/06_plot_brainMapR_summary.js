#!/usr/bin/env node

/*
Create presentation-ready SVG summaries from collected brainMapR outputs.

Inputs:
  summary/matrix_rGM.tsv
  summary/matrix_pvalue_rGM.tsv
  summary/matrix_rGM_CI_lb.tsv
  summary/matrix_rGM_CI_ub.tsv

Outputs:
  summary/figures/figure_1_rGM_heatmap.svg
  summary/figures/figure_2_pvalue_heatmap.svg
  summary/figures/figure_3_top_rGM_forest.svg
  summary/figures/figure_4_batch_qc.svg
  summary/top_associations.tsv
  summary/out_of_range_rGM.tsv

How to run:
  node scripts/06_plot_brainMapR_summary.js --summary-dir /Users/junzhou/Desktop/summary
*/

const fs = require("fs");
const path = require("path");

function parseArgs(argv) {
  const opts = {
    summaryDir: "outputs/batch/summary",
    outDir: null,
  };
  for (let i = 2; i < argv.length; i += 1) {
    const key = argv[i];
    if (key === "--summary-dir") {
      i += 1;
      opts.summaryDir = argv[i];
    } else if (key === "--out-dir") {
      i += 1;
      opts.outDir = argv[i];
    } else if (key === "-h" || key === "--help") {
      console.log("Usage: node scripts/06_plot_brainMapR_summary.js --summary-dir <dir>");
      process.exit(0);
    } else {
      throw new Error(`Unknown argument: ${key}`);
    }
  }
  if (!opts.outDir) {
    opts.outDir = path.join(opts.summaryDir, "figures");
  }
  return opts;
}

function readTsv(file) {
  const text = fs.readFileSync(file, "utf8").trim();
  const lines = text.split(/\r?\n/);
  const header = lines.shift().split("\t");
  return lines.map((line) => {
    const fields = line.split("\t");
    const row = {};
    header.forEach((key, idx) => {
      row[key] = fields[idx] ?? "";
    });
    return row;
  });
}

function matrixToRecords(matrixRows, metricName) {
  const records = [];
  for (const row of matrixRows) {
    const adId = row.ad_id;
    for (const [traitId, value] of Object.entries(row)) {
      if (traitId === "ad_id") continue;
      const numeric = Number(value);
      records.push({
        ad_id: adId,
        trait_id: traitId,
        [metricName]: Number.isFinite(numeric) ? numeric : NaN,
      });
    }
  }
  return records;
}

function inferCategory(traitId) {
  const id = traitId.toLowerCase();
  if (["eid", "sex"].includes(id) || id.endsWith("_pilot")) return "excluded";
  if (/hypert|htn|t2d|diabetes|cholesterol|ldl|bmi|stroke/.test(id)) return "vascular_metabolic";
  if (/depression|mdd|anx|ptsd|bd|scz|phob|psy_/.test(id)) return "psychiatric";
  if (/smoking|pack_year|sleep|insomnia|napping|alcohol|drinker|wine|beer|spirits|pa_/.test(id)) return "lifestyle";
  if (/education|income|imd|social|loneliness|isolation/.test(id)) return "social_socioeconomic";
  if (/frailty|multimorbidity|med_20003/.test(id)) return "frailty_multimorbidity";
  if (/menopause|hrt/.test(id)) return "hormonal_sex_specific";
  if (/infect|lrti|pneumonia|sepsis|ssti|uti|periodontal/.test(id)) return "infection_inflammatory";
  return "other";
}

function bhFdr(records, pKey = "pvalue") {
  const valid = records
    .map((record, idx) => ({ record, idx, p: Number(record[pKey]) }))
    .filter((x) => Number.isFinite(x.p) && x.p >= 0 && x.p <= 1)
    .sort((a, b) => a.p - b.p);
  const m = valid.length;
  let prev = 1;
  for (let i = m - 1; i >= 0; i -= 1) {
    const rank = i + 1;
    const q = Math.min(prev, (valid[i].p * m) / rank);
    valid[i].record.fdr = q;
    prev = q;
  }
  for (const record of records) {
    if (record.fdr === undefined) record.fdr = NaN;
  }
}

function ensureDir(dir) {
  fs.mkdirSync(dir, { recursive: true });
}

function escapeXml(text) {
  return String(text)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function clamp(value, lo, hi) {
  return Math.max(lo, Math.min(hi, value));
}

function lerp(a, b, t) {
  return Math.round(a + (b - a) * t);
}

function divergingColor(value) {
  if (!Number.isFinite(value)) return "#f2f2f2";
  const v = clamp(value, -1, 1);
  if (v < 0) {
    const t = (v + 1);
    return `rgb(${lerp(44, 247, t)},${lerp(123, 247, t)},${lerp(182, 247, t)})`;
  }
  const t = v;
  return `rgb(${lerp(247, 178, t)},${lerp(247, 24, t)},${lerp(247, 43, t)})`;
}

function pColor(p) {
  if (!Number.isFinite(p) || p <= 0) return "#67000d";
  const score = clamp(-Math.log10(p), 0, 20) / 20;
  return `rgb(${lerp(255, 103, score)},${lerp(245, 0, score)},${lerp(235, 31, score)})`;
}

function writeTsv(file, rows, columns) {
  const lines = [columns.join("\t")];
  for (const row of rows) {
    lines.push(columns.map((col) => row[col] ?? "").join("\t"));
  }
  fs.writeFileSync(file, `${lines.join("\n")}\n`);
}

function metricMap(records, key) {
  const map = new Map();
  for (const record of records) {
    map.set(`${record.ad_id}\t${record.trait_id}`, record[key]);
  }
  return map;
}

function categoryLabel(category) {
  const labels = {
    vascular_metabolic: "Vascular/metabolic",
    psychiatric: "Psychiatric",
    lifestyle: "Lifestyle",
    social_socioeconomic: "Social/SES",
    frailty_multimorbidity: "Frailty/multimorbidity",
    infection_inflammatory: "Infection/inflammatory",
    hormonal_sex_specific: "Hormonal/sex-specific",
    other: "Other",
  };
  return labels[category] || category;
}

function orderTraits(traits, records) {
  const maxAbs = new Map();
  for (const trait of traits) maxAbs.set(trait, 0);
  for (const record of records) {
    if (!Number.isFinite(record.rGM)) continue;
    maxAbs.set(record.trait_id, Math.max(maxAbs.get(record.trait_id) || 0, Math.abs(record.rGM)));
  }
  const categoryOrder = [
    "vascular_metabolic",
    "psychiatric",
    "lifestyle",
    "social_socioeconomic",
    "frailty_multimorbidity",
    "infection_inflammatory",
    "hormonal_sex_specific",
    "other",
  ];
  return [...traits].sort((a, b) => {
    const ca = inferCategory(a);
    const cb = inferCategory(b);
    const cDiff = categoryOrder.indexOf(ca) - categoryOrder.indexOf(cb);
    if (cDiff !== 0) return cDiff;
    return (maxAbs.get(b) || 0) - (maxAbs.get(a) || 0) || a.localeCompare(b);
  });
}

function adOrder(rows) {
  const desired = [
    "ADvsHC",
    "MCIvsHC",
    "Conversion1year",
    "Conversion2years",
    "Conversion3years",
    "Conversion4years",
    "Conversion5years",
    "MMSE",
  ];
  return desired.filter((id) => rows.includes(id)).concat(rows.filter((id) => !desired.includes(id)));
}

function heatmapSvg({
  title,
  subtitle,
  adIds,
  traitIds,
  valueMap,
  pMap,
  fdrMap,
  colorFn,
  legendTitle,
  legendTopLabel = "+1 / high",
  legendMiddleLabel = "0",
  legendBottomLabel = "-1 / low",
  valueFormatter,
  outOfRangeMap = new Map(),
}) {
  const cellW = 20;
  const cellH = 28;
  const left = 150;
  const top = 125;
  const right = 230;
  const bottom = 230;
  const width = left + traitIds.length * cellW + right;
  const height = top + adIds.length * cellH + bottom;
  const categoryColors = {
    vascular_metabolic: "#4c78a8",
    psychiatric: "#b279a2",
    lifestyle: "#59a14f",
    social_socioeconomic: "#f28e2b",
    frailty_multimorbidity: "#9c755f",
    infection_inflammatory: "#e15759",
    hormonal_sex_specific: "#76b7b2",
    other: "#bab0ac",
  };
  const chunks = [];
  chunks.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`);
  chunks.push(`<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L6,3 z" fill="#777"/></marker></defs>`);
  chunks.push(`<rect width="100%" height="100%" fill="white"/>`);
  chunks.push(`<style>
    text { font-family: Arial, Helvetica, sans-serif; fill: #222; }
    .title { font-size: 22px; font-weight: 700; }
    .subtitle { font-size: 13px; fill: #555; }
    .axis { font-size: 12px; }
    .small { font-size: 11px; fill: #555; }
    .tiny { font-size: 9px; fill: #555; }
    .categoryLabel { font-size: 7px; fill: #555; }
  </style>`);
  chunks.push(`<text class="title" x="30" y="34">${escapeXml(title)}</text>`);
  chunks.push(`<text class="subtitle" x="30" y="58">${escapeXml(subtitle)}</text>`);

  let categoryStart = 0;
  while (categoryStart < traitIds.length) {
    const category = inferCategory(traitIds[categoryStart]);
    let categoryEnd = categoryStart + 1;
    while (categoryEnd < traitIds.length && inferCategory(traitIds[categoryEnd]) === category) categoryEnd += 1;
    const x = left + categoryStart * cellW;
    const w = (categoryEnd - categoryStart) * cellW;
    chunks.push(`<rect x="${x}" y="${top - 18}" width="${w}" height="8" fill="${categoryColors[category] || "#aaa"}"/>`);
    chunks.push(`<text class="categoryLabel" x="${x + w / 2}" y="${top - 24}" text-anchor="middle">${escapeXml(categoryLabel(category))}</text>`);
    categoryStart = categoryEnd;
  }

  traitIds.forEach((trait, idx) => {
    const x = left + idx * cellW + cellW / 2;
    const y = top + adIds.length * cellH + 18;
    chunks.push(`<text class="axis" transform="translate(${x},${y}) rotate(60)" text-anchor="start">${escapeXml(trait)}</text>`);
  });
  adIds.forEach((ad, rowIdx) => {
    const y = top + rowIdx * cellH + cellH / 2 + 4;
    chunks.push(`<text class="axis" x="${left - 12}" y="${y}" text-anchor="end">${escapeXml(ad)}</text>`);
  });

  adIds.forEach((ad, rowIdx) => {
    traitIds.forEach((trait, colIdx) => {
      const key = `${ad}\t${trait}`;
      const v = valueMap.get(key);
      const p = pMap ? pMap.get(key) : NaN;
      const fdr = fdrMap ? fdrMap.get(key) : NaN;
      const x = left + colIdx * cellW;
      const y = top + rowIdx * cellH;
      chunks.push(`<rect x="${x}" y="${y}" width="${cellW}" height="${cellH}" fill="${colorFn(v, p, fdr)}" stroke="#ffffff" stroke-width="0.6">`);
      chunks.push(`<title>${escapeXml(`${ad} x ${trait}\nvalue=${valueFormatter(v)}\np=${Number.isFinite(p) ? p.toExponential(2) : "NA"}\nFDR=${Number.isFinite(fdr) ? fdr.toExponential(2) : "NA"}`)}</title></rect>`);
      if (Number.isFinite(fdr) && fdr < 0.05) {
        chunks.push(`<circle cx="${x + cellW / 2}" cy="${y + cellH / 2}" r="2.5" fill="#111"/>`);
      }
      if (outOfRangeMap.get(key)) {
        chunks.push(`<text x="${x + cellW / 2}" y="${y + cellH / 2 + 4}" text-anchor="middle" font-size="12" font-weight="700" fill="#111">!</text>`);
      }
    });
  });

  const legendX = left + traitIds.length * cellW + 35;
  const legendY = top;
  chunks.push(`<text class="axis" x="${legendX}" y="${legendY - 18}" font-weight="700">${escapeXml(legendTitle)}</text>`);
  for (let i = 0; i <= 100; i += 1) {
    const t = i / 100;
    const v = -1 + 2 * t;
    chunks.push(`<rect x="${legendX}" y="${legendY + i * 1.8}" width="18" height="2" fill="${colorFn(v, NaN, NaN)}"/>`);
  }
  chunks.push(`<text class="small" x="${legendX + 25}" y="${legendY + 4}">${escapeXml(legendTopLabel)}</text>`);
  chunks.push(`<text class="small" x="${legendX + 25}" y="${legendY + 94}">${escapeXml(legendMiddleLabel)}</text>`);
  chunks.push(`<text class="small" x="${legendX + 25}" y="${legendY + 184}">${escapeXml(legendBottomLabel)}</text>`);
  chunks.push(`<circle cx="${legendX + 5}" cy="${legendY + 230}" r="3" fill="#111"/><text class="small" x="${legendX + 18}" y="${legendY + 234}">FDR &lt; 0.05</text>`);
  chunks.push(`<text x="${legendX}" y="${legendY + 258}" font-size="13" font-weight="700" fill="#111">!</text><text class="small" x="${legendX + 18}" y="${legendY + 258}">|rGM| &gt; 1 (unstable)</text>`);
  chunks.push(`<text class="small" x="30" y="${height - 28}">Cells show collected successful jobs only. Failed/unusable traits are excluded from the matrix outputs.</text>`);
  chunks.push(`</svg>`);
  return chunks.join("\n");
}

function forestSvg(records) {
  const rows = records.slice(0, 30);
  const rowH = 26;
  const left = 390;
  const right = 80;
  const top = 84;
  const bottom = 70;
  const width = 980;
  const height = top + rows.length * rowH + bottom;
  const minX = -1;
  const maxX = 1;
  const scale = (v) => left + ((clamp(v, minX, maxX) - minX) / (maxX - minX)) * (width - left - right);
  const categoryColors = {
    vascular_metabolic: "#4c78a8",
    psychiatric: "#b279a2",
    lifestyle: "#59a14f",
    social_socioeconomic: "#f28e2b",
    frailty_multimorbidity: "#9c755f",
    infection_inflammatory: "#e15759",
    hormonal_sex_specific: "#76b7b2",
    other: "#555",
  };
  const chunks = [];
  chunks.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`);
  chunks.push(`<rect width="100%" height="100%" fill="white"/>`);
  chunks.push(`<style>text{font-family:Arial,Helvetica,sans-serif;fill:#222}.title{font-size:22px;font-weight:700}.sub{font-size:13px;fill:#555}.lab{font-size:12px}.tiny{font-size:10px;fill:#555}</style>`);
  chunks.push(`<text class="title" x="28" y="34">Top stable rGM associations</text>`);
  chunks.push(`<text class="sub" x="28" y="56">Top FDR-significant AD x risk-factor pairs with |rGM| ≤ 1. Error bars show collected 95% CI.</text>`);
  [-1, -0.5, 0, 0.5, 1].forEach((tick) => {
    const x = scale(tick);
    chunks.push(`<line x1="${x}" y1="${top - 18}" x2="${x}" y2="${height - bottom + 6}" stroke="${tick === 0 ? "#222" : "#ddd"}" stroke-width="${tick === 0 ? 1.2 : 1}"/>`);
    chunks.push(`<text class="tiny" x="${x}" y="${height - bottom + 26}" text-anchor="middle">${tick}</text>`);
  });
  rows.forEach((row, idx) => {
    const y = top + idx * rowH;
    const label = `${row.ad_id} x ${row.trait_id}`;
    const color = categoryColors[inferCategory(row.trait_id)] || "#555";
    chunks.push(`<text class="lab" x="${left - 12}" y="${y + 5}" text-anchor="end">${escapeXml(label)}</text>`);
    chunks.push(`<line x1="${scale(row.ci_lb)}" y1="${y}" x2="${scale(row.ci_ub)}" y2="${y}" stroke="${color}" stroke-width="2"/>`);
    chunks.push(`<circle cx="${scale(row.rGM)}" cy="${y}" r="4.5" fill="${color}"><title>${escapeXml(`rGM=${row.rGM.toFixed(3)}\n95% CI ${row.ci_lb.toFixed(3)} to ${row.ci_ub.toFixed(3)}\nFDR=${row.fdr.toExponential(2)}`)}</title></circle>`);
    chunks.push(`<text class="tiny" x="${width - 70}" y="${y + 4}" text-anchor="end">${row.rGM.toFixed(2)}</text>`);
  });
  chunks.push(`<text class="lab" x="${(left + width - right) / 2}" y="${height - 16}" text-anchor="middle">rGM</text>`);
  chunks.push(`</svg>`);
  return chunks.join("\n");
}

function qcSvg() {
  const width = 980;
  const height = 430;
  const boxes = [
    { x: 40, y: 95, w: 190, h: 78, title: "Inputs audited", lines: ["24 AD BWAS files", "71 UKB BWAS files"] },
    { x: 285, y: 95, w: 210, h: 78, title: "Default manifest", lines: ["8 AD maps with confirmed N", "66 usable UKB risk-factor maps"] },
    { x: 550, y: 95, w: 180, h: 78, title: "PBS array", lines: ["528 submitted pairs", "8 x 66 grid"] },
    { x: 785, y: 95, w: 155, h: 78, title: "Completed", lines: ["488 success", "40 failed"] },
    { x: 130, y: 260, w: 300, h: 95, title: "Pre-run exclusions", lines: ["16 AD maps without confirmed N", "2 header-only UKB files", "3 alcohol files corrected in derived inputs"] },
    { x: 555, y: 260, w: 300, h: 95, title: "Trait-level failures", lines: ["5 UKB traits x 8 AD maps", "b=0, se=0, p=NA", "excluded unless regenerated upstream"] },
  ];
  const chunks = [];
  chunks.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`);
  chunks.push(`<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L6,3 z" fill="#777"/></marker></defs>`);
  chunks.push(`<rect width="100%" height="100%" fill="white"/>`);
  chunks.push(`<style>text{font-family:Arial,Helvetica,sans-serif;fill:#222}.title{font-size:22px;font-weight:700}.sub{font-size:13px;fill:#555}.boxTitle{font-size:15px;font-weight:700}.line{font-size:13px;fill:#333}</style>`);
  chunks.push(`<text class="title" x="36" y="38">brainMapR batch QC summary</text>`);
  chunks.push(`<text class="sub" x="36" y="62">Pipeline status after full AD x risk-factor PBS array.</text>`);
  boxes.forEach((box) => {
    chunks.push(`<rect x="${box.x}" y="${box.y}" width="${box.w}" height="${box.h}" rx="6" fill="#f8f9fb" stroke="#c7cbd1"/>`);
    chunks.push(`<text class="boxTitle" x="${box.x + 14}" y="${box.y + 26}">${escapeXml(box.title)}</text>`);
    box.lines.forEach((line, idx) => {
      chunks.push(`<text class="line" x="${box.x + 14}" y="${box.y + 50 + idx * 18}">${escapeXml(line)}</text>`);
    });
  });
  boxes.slice(0, 3).forEach((box, idx) => {
    const next = boxes[idx + 1];
    const y = box.y + box.h / 2;
    chunks.push(`<line x1="${box.x + box.w}" y1="${y}" x2="${next.x}" y2="${y}" stroke="#777" stroke-width="1.5" marker-end="url(#arrow)"/>`);
  });
  chunks.push(`</svg>`);
  return chunks.join("\n");
}

function main() {
  const opts = parseArgs(process.argv);
  ensureDir(opts.outDir);

  const rRows = readTsv(path.join(opts.summaryDir, "matrix_rGM.tsv"));
  const pRows = readTsv(path.join(opts.summaryDir, "matrix_pvalue_rGM.tsv"));
  const lbRows = readTsv(path.join(opts.summaryDir, "matrix_rGM_CI_lb.tsv"));
  const ubRows = readTsv(path.join(opts.summaryDir, "matrix_rGM_CI_ub.tsv"));

  const rRecords = matrixToRecords(rRows, "rGM");
  const pRecords = matrixToRecords(pRows, "pvalue");
  const lbRecords = matrixToRecords(lbRows, "ci_lb");
  const ubRecords = matrixToRecords(ubRows, "ci_ub");
  const pByKey = metricMap(pRecords, "pvalue");
  const lbByKey = metricMap(lbRecords, "ci_lb");
  const ubByKey = metricMap(ubRecords, "ci_ub");

  const records = rRecords.map((record) => {
    const key = `${record.ad_id}\t${record.trait_id}`;
    return {
      ...record,
      pvalue: pByKey.get(key),
      ci_lb: lbByKey.get(key),
      ci_ub: ubByKey.get(key),
      category: inferCategory(record.trait_id),
    };
  });
  bhFdr(records, "pvalue");

  const adIds = adOrder([...new Set(records.map((x) => x.ad_id))]);
  const traitIds = orderTraits([...new Set(records.map((x) => x.trait_id))], records);
  const rMap = metricMap(records, "rGM");
  const pMap = metricMap(records, "pvalue");
  const fdrMap = metricMap(records, "fdr");
  const oorMap = new Map();
  records.forEach((record) => {
    if (Number.isFinite(record.rGM) && Math.abs(record.rGM) > 1) {
      oorMap.set(`${record.ad_id}\t${record.trait_id}`, true);
    }
  });

  const heatmap = heatmapSvg({
    title: "AD x risk-factor grey-matter correlations (rGM)",
    subtitle: "Color shows rGM clipped to [-1, 1]; dot marks FDR < 0.05; ! marks out-of-range estimates kept in the TSV.",
    adIds,
    traitIds,
    valueMap: rMap,
    pMap,
    fdrMap,
    colorFn: (v) => divergingColor(v),
    legendTitle: "rGM",
    legendTopLabel: "+1",
    legendMiddleLabel: "0",
    legendBottomLabel: "-1",
    valueFormatter: (v) => (Number.isFinite(v) ? v.toFixed(3) : "NA"),
    outOfRangeMap: oorMap,
  });
  fs.writeFileSync(path.join(opts.outDir, "figure_1_rGM_heatmap.svg"), heatmap);

  const pHeatmap = heatmapSvg({
    title: "AD x risk-factor rGM significance",
    subtitle: "Color shows -log10(p) clipped at 20; dot marks FDR < 0.05.",
    adIds,
    traitIds,
    valueMap: pMap,
    pMap,
    fdrMap,
    colorFn: (v, p) => {
      if (Number.isFinite(p)) return pColor(p);
      const negLog10P = ((v + 1) / 2) * 20;
      return pColor(10 ** -negLog10P);
    },
    legendTitle: "-log10(p)",
    legendTopLabel: "20",
    legendMiddleLabel: "10",
    legendBottomLabel: "0",
    valueFormatter: (v) => (Number.isFinite(v) ? v.toExponential(2) : "NA"),
  });
  fs.writeFileSync(path.join(opts.outDir, "figure_2_pvalue_heatmap.svg"), pHeatmap);

  const top = records
    .filter((x) =>
      Number.isFinite(x.rGM) &&
      Number.isFinite(x.pvalue) &&
      Number.isFinite(x.fdr) &&
      Number.isFinite(x.ci_lb) &&
      Number.isFinite(x.ci_ub) &&
      Math.abs(x.rGM) <= 1 &&
      x.fdr < 0.05
    )
    .sort((a, b) => Math.abs(b.rGM) - Math.abs(a.rGM));
  fs.writeFileSync(path.join(opts.outDir, "figure_3_top_rGM_forest.svg"), forestSvg(top));

  fs.writeFileSync(path.join(opts.outDir, "figure_4_batch_qc.svg"), qcSvg());

  const topRows = top.slice(0, 50).map((x) => ({
    ad_id: x.ad_id,
    trait_id: x.trait_id,
    category: x.category,
    rGM: x.rGM,
    ci_lb: x.ci_lb,
    ci_ub: x.ci_ub,
    pvalue: x.pvalue,
    fdr: x.fdr,
  }));
  writeTsv(path.join(opts.summaryDir, "top_associations.tsv"), topRows, [
    "ad_id",
    "trait_id",
    "category",
    "rGM",
    "ci_lb",
    "ci_ub",
    "pvalue",
    "fdr",
  ]);

  const outOfRange = records
    .filter((x) => Number.isFinite(x.rGM) && Math.abs(x.rGM) > 1)
    .sort((a, b) => Math.abs(b.rGM) - Math.abs(a.rGM));
  writeTsv(path.join(opts.summaryDir, "out_of_range_rGM.tsv"), outOfRange, [
    "ad_id",
    "trait_id",
    "category",
    "rGM",
    "ci_lb",
    "ci_ub",
    "pvalue",
    "fdr",
  ]);

  const report = [
    `Figures written to: ${opts.outDir}`,
    `AD maps: ${adIds.length}`,
    `Risk-factor maps: ${traitIds.length}`,
    `Collected pairs: ${records.length}`,
    `FDR-significant pairs: ${records.filter((x) => Number.isFinite(x.fdr) && x.fdr < 0.05).length}`,
    `Out-of-range rGM estimates: ${outOfRange.length}`,
    `Top stable associations written: ${topRows.length}`,
  ].join("\n");
  fs.writeFileSync(path.join(opts.summaryDir, "figure_generation_report.txt"), `${report}\n`);
  console.log(report);
}

main();
