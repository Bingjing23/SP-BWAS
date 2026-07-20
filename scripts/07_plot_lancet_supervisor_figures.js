#!/usr/bin/env node

/*
Create clean supervisor-facing figures comparing Lancet 2024 dementia
risk-factor priority with current brainMapR map-level results.

Inputs:
  Values are defined in the `riskFactors` table below, using:
    - Lancet 2024 rounded weighted PAF values
    - Current AVERAGE reference-panel brainMapR max abs(rGM) summaries

Outputs:
  outputs/lancet_2024_supervisor_figures_clean/
    figure_1_lancet_rank_lollipop.svg
    figure_2_lancet_paf_rgm_quadrant.svg
    README.md

If `rsvg-convert` is available, PNG versions are also written.

How to run:
  node scripts/07_plot_lancet_supervisor_figures.js
*/

const fs = require("fs");
const path = require("path");
const childProcess = require("child_process");

const projectRoot = path.resolve(__dirname, "..");
const outputDir = path.join(projectRoot, "outputs", "lancet_2024_supervisor_figures_clean");

const weakThreshold = 0.2;

const riskFactors = [
  { rank: 1, factor: "Hearing loss", paf: 7, trait: "Not tested", maxAbsRgm: null, status: "not_tested" },
  { rank: 2, factor: "High LDL cholesterol", paf: 7, trait: "ldl_direct", maxAbsRgm: 0.444, status: "clear" },
  { rank: 3, factor: "Less education", paf: 5, trait: "age_completed_education", maxAbsRgm: 0.176, status: "weak" },
  { rank: 4, factor: "Social isolation", paf: 5, trait: "loneliness_isolation", maxAbsRgm: 0.700, status: "clear" },
  { rank: 5, factor: "Depression", paf: 3, trait: "majorDepression", maxAbsRgm: 0.629, status: "clear" },
  { rank: 6, factor: "Traumatic brain injury", paf: 3, trait: "Not tested", maxAbsRgm: null, status: "not_tested" },
  { rank: 7, factor: "Air pollution", paf: 3, trait: "Not tested", maxAbsRgm: null, status: "not_tested" },
  { rank: 8, factor: "Physical inactivity", paf: 2, trait: "pa_* proxies", maxAbsRgm: 0.550, status: "mixed" },
  { rank: 9, factor: "Smoking", paf: 2, trait: "smoking_ever / pack_years", maxAbsRgm: 0.207, status: "modest" },
  { rank: 10, factor: "Diabetes", paf: 2, trait: "T2D", maxAbsRgm: 0.067, status: "weak" },
  { rank: 11, factor: "Hypertension", paf: 2, trait: "hyperTension", maxAbsRgm: 0.133, status: "weak" },
  { rank: 12, factor: "Untreated vision loss", paf: 2, trait: "Not tested", maxAbsRgm: null, status: "not_tested" },
  { rank: 13, factor: "Obesity", paf: 1, trait: "bmi", maxAbsRgm: 0.435, status: "clear" },
  { rank: 14, factor: "Excessive alcohol", paf: 1, trait: "alcohol proxies", maxAbsRgm: 0.166, status: "weak" },
];

const palette = {
  clear: "#15803D",
  weak: "#C2413A",
  mixed: "#7C4D9E",
  modest: "#C77700",
  not_tested: "#8B929C",
  text: "#1F2937",
  muted: "#5E6A78",
  grid: "#D7DDE5",
  panel: "#F7F8FA",
  weakFill: "#FDECEC",
  notTestedFill: "#F1F3F5",
  pafBar: "#9DB7E8",
  axis: "#334155",
};

function esc(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function text(x, y, value, options = {}) {
  const size = options.size || 13;
  const weight = options.weight || 400;
  const fill = options.fill || palette.text;
  const anchor = options.anchor || "start";
  return [
    `<text x="${x}" y="${y}"`,
    `font-family="Arial, Helvetica, sans-serif"`,
    `font-size="${size}" font-weight="${weight}" fill="${fill}"`,
    `text-anchor="${anchor}">${esc(value)}</text>`,
  ].join(" ");
}

function rect(x, y, width, height, fill, stroke = "none", rx = 0) {
  return `<rect x="${x}" y="${y}" width="${width}" height="${height}" fill="${fill}" stroke="${stroke}" rx="${rx}"/>`;
}

function line(x1, y1, x2, y2, stroke = palette.grid, width = 1, dash = "") {
  return `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${stroke}" stroke-width="${width}"${dash ? ` stroke-dasharray="${dash}"` : ""}/>`;
}

function circle(cx, cy, radius, fill, stroke = "white", strokeWidth = 1.5) {
  return `<circle cx="${cx}" cy="${cy}" r="${radius}" fill="${fill}" stroke="${stroke}" stroke-width="${strokeWidth}"/>`;
}

function statusLabel(status) {
  return {
    clear: "Clearer map signal",
    weak: "Weak in our maps",
    mixed: "Mixed / QC caution",
    modest: "Modest",
    not_tested: "Not tested",
  }[status];
}

function valueLabel(value) {
  return value === null ? "NA" : value.toFixed(3);
}

function svgDocument(width, height, content) {
  return [
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`,
    rect(0, 0, width, height, "#FFFFFF"),
    ...content,
    "</svg>",
  ].join("\n");
}

function makeFigure1() {
  const width = 1600;
  const height = 900;
  const top = 138;
  const rowHeight = 42;
  const rowLeft = 36;
  const rowRight = width - 36;

  const riskX = 70;
  const traitX = 385;
  const pafX = 700;
  const pafWidth = 260;
  const rgmX = 1060;
  const rgmWidth = 300;
  const statusX = 1430;

  const pafScale = (value) => pafX + (value / 7) * pafWidth;
  const rgmScale = (value) => rgmX + (value / 0.75) * rgmWidth;

  const out = [];
  out.push(text(42, 42, "Lancet 2024 risk-factor priority vs current brainMapR signal", { size: 25, weight: 700 }));
  out.push(text(42, 70, "Ordered by Lancet weighted PAF. Red rows indicate abs(rGM) < 0.20 in our current map-level results.", { size: 15, fill: palette.muted }));

  out.push(text(riskX, 112, "Lancet risk factor", { size: 14, weight: 700 }));
  out.push(text(traitX, 112, "Closest current trait", { size: 14, weight: 700 }));
  out.push(text(pafX, 112, "Weighted PAF (%)", { size: 14, weight: 700 }));
  out.push(text(rgmX, 112, "Max abs(rGM)", { size: 14, weight: 700 }));
  out.push(text(statusX, 112, "Status", { size: 14, weight: 700 }));

  [0, 2, 4, 6].forEach((tick) => {
    const x = pafScale(tick);
    out.push(line(x, top - 18, x, top + riskFactors.length * rowHeight - 8));
    out.push(text(x, top + riskFactors.length * rowHeight + 20, tick, { size: 11, fill: palette.muted, anchor: "middle" }));
  });

  [0, 0.2, 0.4, 0.6].forEach((tick) => {
    const x = rgmScale(tick);
    out.push(line(x, top - 18, x, top + riskFactors.length * rowHeight - 8));
    out.push(text(x, top + riskFactors.length * rowHeight + 20, tick.toFixed(1), { size: 11, fill: palette.muted, anchor: "middle" }));
  });

  out.push(line(rgmScale(weakThreshold), top - 20, rgmScale(weakThreshold), top + riskFactors.length * rowHeight - 8, palette.weak, 2, "5 5"));
  out.push(text(rgmScale(weakThreshold) + 8, 126, "weak threshold", { size: 11, fill: palette.weak, weight: 700 }));

  riskFactors.forEach((row, index) => {
    const y = top + index * rowHeight;
    if (row.status === "weak") {
      out.push(rect(rowLeft, y - 24, rowRight - rowLeft, rowHeight - 4, palette.weakFill, "none", 7));
    } else if (row.status === "not_tested") {
      out.push(rect(rowLeft, y - 24, rowRight - rowLeft, rowHeight - 4, palette.notTestedFill, "none", 7));
    }

    out.push(text(riskX, y, `${row.rank}. ${row.factor}`, { size: 14, weight: row.status === "weak" ? 700 : 500 }));
    out.push(text(traitX, y, row.trait, { size: 12, fill: palette.muted }));

    out.push(rect(pafX, y - 13, pafScale(row.paf) - pafX, 16, palette.pafBar, "none", 2));
    out.push(text(pafScale(row.paf) + 8, y, `${row.paf}%`, { size: 12, fill: palette.text }));

    if (row.maxAbsRgm === null) {
      out.push(text(rgmX, y, "not available", { size: 12, fill: palette.not_tested }));
    } else {
      out.push(line(rgmX, y - 5, rgmScale(row.maxAbsRgm), y - 5, palette[row.status], 5));
      out.push(circle(rgmScale(row.maxAbsRgm), y - 5, 7, palette[row.status]));
      out.push(text(rgmScale(row.maxAbsRgm) + 10, y, valueLabel(row.maxAbsRgm), {
        size: 12,
        fill: palette[row.status],
        weight: 700,
      }));
    }

    out.push(text(statusX, y, statusLabel(row.status), {
      size: 12,
      fill: palette[row.status],
      weight: row.status === "weak" ? 700 : 500,
    }));
  });

  const legendY = height - 42;
  let legendX = 50;
  [
    ["Weak in our maps", palette.weak],
    ["Clearer signal", palette.clear],
    ["Modest / mixed", palette.mixed],
    ["Not tested", palette.not_tested],
  ].forEach(([label, color]) => {
    out.push(circle(legendX, legendY - 4, 6, color));
    out.push(text(legendX + 14, legendY, label, { size: 12, fill: palette.muted }));
    legendX += 180;
  });

  out.push(text(width - 46, height - 18, "Metric: rGM from brainMapR sumR2_regression_bivariate; reference panel: AVERAGE", {
    size: 11,
    fill: palette.muted,
    anchor: "end",
  }));

  return svgDocument(width, height, out);
}

function makeFigure2() {
  const width = 1600;
  const height = 940;
  const margin = { left: 92, right: 455, top: 122, bottom: 92 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const xMax = 7.2;
  const yMax = 0.75;
  const xScale = (value) => margin.left + (value / xMax) * plotWidth;
  const yScale = (value) => margin.top + plotHeight - (value / yMax) * plotHeight;

  const out = [];
  out.push(text(44, 44, "Lancet-prioritised factors with weak current map-level signal", { size: 25, weight: 700 }));
  out.push(text(44, 72, "Numbers identify risk factors. The red zone marks abs(rGM) < 0.20; not-tested factors are placed near the baseline.", { size: 15, fill: palette.muted }));

  out.push(rect(xScale(0), yScale(weakThreshold), xScale(xMax) - xScale(0), yScale(0) - yScale(weakThreshold), palette.weakFill, "#F2B8B5", 3));
  out.push(text(xScale(3.6), yScale(0.08), "Weak map-level signal zone", { size: 14, fill: palette.weak, weight: 700, anchor: "middle" }));

  [0, 1, 2, 3, 4, 5, 6, 7].forEach((tick) => {
    const x = xScale(tick);
    out.push(line(x, margin.top, x, margin.top + plotHeight));
    out.push(text(x, margin.top + plotHeight + 24, tick, { size: 12, fill: palette.muted, anchor: "middle" }));
  });

  [0, 0.2, 0.4, 0.6].forEach((tick) => {
    const y = yScale(tick);
    out.push(line(margin.left, y, margin.left + plotWidth, y));
    out.push(text(margin.left - 12, y + 4, tick.toFixed(1), { size: 12, fill: palette.muted, anchor: "end" }));
  });

  out.push(line(margin.left, margin.top + plotHeight, margin.left + plotWidth, margin.top + plotHeight, palette.axis, 1.5));
  out.push(line(margin.left, margin.top, margin.left, margin.top + plotHeight, palette.axis, 1.5));
  out.push(text(margin.left + plotWidth / 2, height - 38, "Lancet 2024 weighted PAF (%)", { size: 15, weight: 700, anchor: "middle" }));
  out.push(`<text x="26" y="${margin.top + plotHeight / 2}" font-family="Arial, Helvetica, sans-serif" font-size="15" font-weight="700" fill="${palette.text}" transform="rotate(-90 26 ${margin.top + plotHeight / 2})" text-anchor="middle">Our maximum abs(rGM) across AD-related maps</text>`);

  riskFactors.forEach((row) => {
    const yValue = row.maxAbsRgm === null ? 0.025 : Math.min(row.maxAbsRgm, yMax - 0.02);
    const x = xScale(row.paf);
    const y = yScale(yValue);
    const radius = row.status === "weak" ? 12 : 11;

    if (row.maxAbsRgm === null) {
      out.push(circle(x, y, radius, "#FFFFFF", palette.not_tested, 2));
      out.push(text(x, y + 4, row.rank, { size: 10, fill: palette.not_tested, weight: 700, anchor: "middle" }));
    } else {
      out.push(circle(x, y, radius, palette[row.status]));
      out.push(text(x, y + 4, row.rank, { size: 10, fill: "#FFFFFF", weight: 700, anchor: "middle" }));
    }
  });

  const tableX = 1190;
  let tableY = 130;
  out.push(text(tableX, tableY - 24, "Risk-factor key", { size: 17, weight: 700 }));
  out.push(text(tableX, tableY, "No.", { size: 11, weight: 700, fill: palette.muted }));
  out.push(text(tableX + 42, tableY, "Factor", { size: 11, weight: 700, fill: palette.muted }));
  out.push(text(tableX + 265, tableY, "Max abs(rGM)", { size: 11, weight: 700, fill: palette.muted }));
  tableY += 18;

  riskFactors.forEach((row) => {
    const fill = palette[row.status];
    out.push(circle(tableX + 8, tableY - 4, 7, row.maxAbsRgm === null ? "#FFFFFF" : fill, row.maxAbsRgm === null ? fill : "#FFFFFF", 1.5));
    out.push(text(tableX + 8, tableY, row.rank, { size: 8, fill: row.maxAbsRgm === null ? fill : "#FFFFFF", weight: 700, anchor: "middle" }));
    out.push(text(tableX + 42, tableY, row.factor, { size: 11, fill: row.status === "weak" ? palette.weak : palette.text, weight: row.status === "weak" ? 700 : 400 }));
    out.push(text(tableX + 285, tableY, valueLabel(row.maxAbsRgm), { size: 11, fill: row.status === "weak" ? palette.weak : palette.muted, weight: row.status === "weak" ? 700 : 400 }));
    tableY += 22;
  });

  tableY += 20;
  out.push(text(tableX, tableY, "Main weak mismatches", { size: 14, weight: 700 }));
  tableY += 24;
  ["Diabetes", "Hypertension", "Less education", "Excessive alcohol"].forEach((label) => {
    out.push(text(tableX, tableY, `- ${label}`, { size: 12, fill: palette.weak }));
    tableY += 22;
  });

  tableY += 14;
  out.push(text(tableX, tableY, "Not directly tested", { size: 14, weight: 700 }));
  tableY += 24;
  ["Hearing loss", "TBI", "Air pollution", "Vision loss"].forEach((label) => {
    out.push(text(tableX, tableY, `- ${label}`, { size: 12, fill: palette.not_tested }));
    tableY += 22;
  });

  out.push(text(width - 42, height - 18, "Interpretation is map-level similarity only; no causal interpretation.", {
    size: 11,
    fill: palette.muted,
    anchor: "end",
  }));

  return svgDocument(width, height, out);
}

function writeOutput() {
  fs.mkdirSync(outputDir, { recursive: true });
  const figure1Svg = path.join(outputDir, "figure_1_lancet_rank_lollipop.svg");
  const figure2Svg = path.join(outputDir, "figure_2_lancet_paf_rgm_quadrant.svg");
  fs.writeFileSync(figure1Svg, makeFigure1());
  fs.writeFileSync(figure2Svg, makeFigure2());

  const readme = [
    "# Clean Lancet 2024 Supervisor Figures",
    "",
    "Generated by `node scripts/07_plot_lancet_supervisor_figures.js`.",
    "",
    "Files:",
    "- `figure_1_lancet_rank_lollipop.svg`",
    "- `figure_2_lancet_paf_rgm_quadrant.svg`",
    "",
    "If `rsvg-convert` is available, PNG versions are generated automatically.",
    "",
    "Weak map-level signal is defined as `abs(rGM) < 0.20`.",
    "",
  ].join("\n");
  fs.writeFileSync(path.join(outputDir, "README.md"), readme);

  const hasRsvg = childProcess.spawnSync("which", ["rsvg-convert"], { encoding: "utf8" });
  if (hasRsvg.status === 0) {
    childProcess.spawnSync("rsvg-convert", ["-w", "1800", "-o", path.join(outputDir, "figure_1_lancet_rank_lollipop.png"), figure1Svg], { stdio: "inherit" });
    childProcess.spawnSync("rsvg-convert", ["-w", "1800", "-o", path.join(outputDir, "figure_2_lancet_paf_rgm_quadrant.png"), figure2Svg], { stdio: "inherit" });
  }

  console.log(`Wrote figures to: ${outputDir}`);
}

writeOutput();
