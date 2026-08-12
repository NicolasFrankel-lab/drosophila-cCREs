"""Integrated, deterministic QC reporting for the regulatory catalog."""

from __future__ import annotations

import base64
import csv
import gzip
import html
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any

from .activity import sha256_file, write_json_if_changed
from .activity_tmm import _atomic_text_if_changed


REPORT_SOURCE_PATTERNS = (
    ("manifest", "provenance/manifests/*.tsv"),
    ("resolved_config", "provenance/resolved_config.json"),
    ("qc_metrics", "qc/metrics.json"),
    ("master_metrics", "atac/master/master_dhs.json"),
    ("tss_profile", "qc/tss/*.profile.png"),
    ("fragment_histogram", "qc/fragments/*.histogram.tsv"),
    ("cross_correlation", "qc/chip/*.cross_correlation.txt"),
    ("multiqc_report", "qc/multiqc/multiqc_report.html"),
)


def report_source_record(path: Path, *, kind: str, source_root: Path) -> dict[str, Any]:
    """Return one checksummed report input record."""

    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Report source file does not exist: {path}")
    return {
        "path": path,
        "sha256": sha256_file(path),
        "kind": kind,
        "source_root": source_root.resolve(),
    }


def discover_report_source_files(source_roots: list[Path]) -> list[dict[str, Any]]:
    """Discover the bounded set of upstream artifacts used by the report."""

    records: dict[Path, dict[str, Any]] = {}
    for source_root in sorted({path.resolve() for path in source_roots}):
        if not source_root.is_dir():
            raise FileNotFoundError(
                f"Report source root does not exist or is not a directory: {source_root}"
            )
        for kind, pattern in REPORT_SOURCE_PATTERNS:
            for path in sorted(source_root.glob(pattern)):
                if path.is_file():
                    records[path.resolve()] = report_source_record(
                        path,
                        kind=kind,
                        source_root=source_root,
                    )
    return [records[path] for path in sorted(records)]


def _read_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _format_value(value: Any) -> str:
    if value is None or value == "":
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "—"
        if abs(value) >= 1000:
            return f"{value:,.0f}"
        return f"{value:.4g}"
    return str(value)


def _table(headers: list[str], rows: list[list[Any]], *, compact: bool = False) -> str:
    class_name = "compact" if compact else ""
    head = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>"
        + "".join(f"<td>{html.escape(_format_value(value))}</td>" for value in row)
        + "</tr>"
        for row in rows
    )
    if not rows:
        body = f'<tr><td colspan="{len(headers)}" class="muted">No records</td></tr>'
    return f'<table class="{class_name}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def _cards(values: list[tuple[str, Any]]) -> str:
    return '<div class="cards">' + "".join(
        '<div class="card">'
        f'<div class="card-value">{html.escape(_format_value(value))}</div>'
        f'<div class="card-label">{html.escape(label)}</div>'
        "</div>"
        for label, value in values
    ) + "</div>"


def _bar_chart(
    title: str,
    values: list[tuple[str, float]],
    *,
    value_format: str = "number",
    color: str = "#3976a8",
) -> str:
    if not values:
        return '<p class="muted">No values available.</p>'
    panel_width = 760
    row_height = 24
    left = 145
    right = 70
    top = 38
    height = top + row_height * len(values) + 28
    plot_width = panel_width - left - right
    maximum = max(value for _label, value in values)
    if maximum <= 0:
        maximum = 1.0
    parts = [
        f'<svg class="chart" viewBox="0 0 {panel_width} {height}" role="img" aria-label="{html.escape(title)}">',
        f'<text x="{panel_width / 2}" y="19" text-anchor="middle" class="chart-title">{html.escape(title)}</text>',
    ]
    for index, (label, value) in enumerate(values):
        y = top + index * row_height
        width = max(0.0, value) / maximum * plot_width
        displayed = (
            f"{value:.3f}" if value_format == "fraction" else _format_value(value)
        )
        parts.extend(
            [
                f'<text x="{left - 7}" y="{y + 13}" text-anchor="end">{html.escape(label)}</text>',
                f'<rect x="{left}" y="{y}" width="{plot_width}" height="16" fill="#edf1f5"/>',
                f'<rect x="{left}" y="{y}" width="{width:.2f}" height="16" rx="2" fill="{color}"/>',
                f'<text x="{left + width + 5:.2f}" y="{y + 13}">{html.escape(displayed)}</text>',
            ]
        )
    parts.append("</svg>")
    return "".join(parts)


def _stacked_chart(
    title: str,
    contexts: list[str],
    counts: dict[str, dict[str, int]],
) -> str:
    categories = (
        "promoter_associated",
        "proximal_enhancer_like",
        "distal_enhancer_like",
        "unclassified_no_tss_on_contig",
    )
    colors = ("#7851a9", "#df8f44", "#3b8f6f", "#8b98a8")
    totals = {context: sum(counts.get(context, {}).values()) for context in contexts}
    maximum = max(totals.values(), default=1) or 1
    width, left, right, top, row_height = 780, 90, 70, 45, 29
    plot_width = width - left - right
    height = top + len(contexts) * row_height + 52
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(title)}">',
        f'<text x="{width / 2}" y="19" text-anchor="middle" class="chart-title">{html.escape(title)}</text>',
    ]
    legend_x = left
    for category, color in zip(categories, colors):
        label = category.replace("_", " ")
        parts.extend(
            [
                f'<rect x="{legend_x}" y="27" width="9" height="9" fill="{color}"/>',
                f'<text x="{legend_x + 13}" y="35">{html.escape(label)}</text>',
            ]
        )
        legend_x += 150
    for index, context in enumerate(contexts):
        y = top + index * row_height
        parts.append(f'<text x="{left - 7}" y="{y + 15}" text-anchor="end">{html.escape(context)}</text>')
        x = left
        for category, color in zip(categories, colors):
            value = counts.get(context, {}).get(category, 0)
            bar_width = value / maximum * plot_width
            parts.append(
                f'<rect x="{x:.2f}" y="{y}" width="{bar_width:.2f}" height="18" fill="{color}"/>'
            )
            x += bar_width
        parts.append(f'<text x="{x + 5:.2f}" y="{y + 15}">{totals[context]:,}</text>')
    parts.append("</svg>")
    return "".join(parts)


def _data_uri(path: Path) -> str:
    media_type = {
        ".png": "image/png",
        ".svg": "image/svg+xml",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
    }.get(path.suffix.lower(), "application/octet-stream")
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{media_type};base64,{encoded}"


def _source_records(source_files: list[dict[str, Any]]) -> list[dict[str, Any]]:
    records = []
    seen = set()
    for source in source_files:
        path = Path(str(source["path"])).resolve()
        if path in seen:
            raise ValueError(f"Duplicate report source path: {path}")
        seen.add(path)
        if not path.is_file():
            raise FileNotFoundError(f"Report source file does not exist: {path}")
        observed = sha256_file(path)
        if observed != source["sha256"]:
            raise ValueError(f"Report source checksum differs: {path}")
        records.append({**source, "path": path, "observed_sha256": observed})
    return records


def _qc_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for record in records:
        if record["kind"] != "qc_metrics":
            continue
        payload = _read_json(record["path"])
        if not isinstance(payload, list):
            raise ValueError(f"QC metrics must be a list: {record['path']}")
        for row in payload:
            sample = str(row.get("sample", ""))
            if not sample:
                raise ValueError(f"QC metric lacks sample: {record['path']}")
            if sample in rows and rows[sample] != row:
                raise ValueError(f"Conflicting upstream QC metrics for {sample}")
            rows[sample] = row
    return [rows[sample] for sample in sorted(rows)]


def _cross_correlation_rows(records: list[dict[str, Any]]) -> list[list[Any]]:
    rows = []
    for record in records:
        if record["kind"] != "cross_correlation":
            continue
        lines = [
            line.split("\t")
            for line in record["path"].read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        for fields in lines:
            if len(fields) < 11:
                raise ValueError(
                    f"Unexpected phantompeakqualtools row: {record['path']}"
                )
            sample = Path(fields[0]).name.removesuffix(".final.bam")
            rows.append(
                [sample, int(fields[1]), fields[2], float(fields[-3]), float(fields[-2]), int(fields[-1])]
            )
    return sorted(rows, key=lambda row: row[0])


def _default_pdf_renderer(html_text: str, output_path: Path) -> str:
    try:
        import weasyprint
    except ImportError as error:  # pragma: no cover - exercised in rule environment
        raise RuntimeError("WeasyPrint is required for the final PDF report") from error
    os.environ.setdefault("SOURCE_DATE_EPOCH", "0")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{output_path.name}.",
            suffix=".pdf",
            dir=output_path.parent,
            delete=False,
        ) as handle:
            temporary_name = handle.name
        weasyprint.HTML(string=html_text, base_url=str(Path.cwd())).write_pdf(
            temporary_name
        )
        temporary = Path(temporary_name)
        if output_path.is_file() and output_path.read_bytes() == temporary.read_bytes():
            temporary.unlink()
        else:
            os.replace(temporary, output_path)
        temporary_name = None
    finally:
        if temporary_name and Path(temporary_name).exists():
            Path(temporary_name).unlink()
    return str(weasyprint.__version__)


def build_integrated_qc_report(
    *,
    config: dict[str, Any],
    source_files: list[dict[str, Any]],
    current_files: dict[str, Path],
    output_html: Path,
    output_pdf: Path,
    output_metrics: Path,
    pdf_renderer=None,
) -> dict[str, Any]:
    """Build an integrated HTML/PDF report from frozen upstream artifacts."""

    if config.get("assay") != "activity" or not isinstance(config.get("activity"), dict):
        raise ValueError("Integrated report requires an activity workflow configuration")
    records = _source_records(source_files)
    for name, path in current_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"Current report artifact {name} does not exist: {path}")

    activity = config["activity"]
    contexts = [str(context) for context in activity["contexts"]]
    library_metadata = {
        str(library["id"]): library for library in activity["libraries"]
    }
    excluded = list(config.get("provenance", {}).get("excluded_activity_libraries", []))
    for library in excluded:
        library_metadata[str(library["id"])] = library

    master_metrics = _read_json(Path(activity["master"]["stats_json"]))
    activity_metrics = _read_json(current_files["activity_metrics"])
    catalog_metrics = _read_json(current_files["regulatory_element_metrics"])
    factors = _read_tsv(current_files["normalization_factors"])
    summary = _read_tsv(current_files["regulatory_element_summary"])
    nearest_tss_metrics = _read_json(current_files["nearest_tss_enhancer_metrics"])
    contact_metrics = (
        _read_json(current_files["contact_graph_metrics"])
        if "contact_graph_metrics" in current_files
        else None
    )
    qc_rows = _qc_rows(records)
    cross_correlation = _cross_correlation_rows(records)

    warnings = []
    for library in excluded:
        warnings.append(f"{library['id']} rejected: {library.get('reason', 'no reason recorded')}")
    accepted_ids = {str(library["id"]) for library in activity["libraries"]}
    qc_ids = {str(row["sample"]) for row in qc_rows}
    missing_qc = sorted(accepted_ids - qc_ids)
    if missing_qc:
        warnings.append("No upstream qc/metrics.json row for: " + ", ".join(missing_qc))

    qc_table_rows = []
    for row in qc_rows:
        sample = str(row["sample"])
        metadata = library_metadata.get(sample, {})
        qc_table_rows.append(
            [
                sample,
                metadata.get("context", ""),
                "H3K27ac" if row.get("assay") == "chip_histone" else str(row.get("assay", "")).upper(),
                row.get("layout"),
                row.get("mapped_alignments"),
                row.get("frip"),
                row.get("peak_count"),
                metadata.get("qc_status", "upstream QC"),
            ]
        )

    active_counts: dict[str, dict[str, int]] = {context: {} for context in contexts}
    for row in summary:
        if row["h3k27ac_active"] != "1":
            continue
        context = row["context"]
        element_class = row["regulatory_class"]
        active_counts.setdefault(context, {})[element_class] = (
            active_counts.setdefault(context, {}).get(element_class, 0)
            + int(row["element_n"])
        )

    factor_rows = [
        [
            row["library_id"],
            row["context"],
            row["assay"],
            int(row["total_units"]),
            float(row["tmm_normalization_factor"]),
            float(row["effective_library_size"]),
            int(row["feature_count"]),
        ]
        for row in factors
    ]
    provenance = config.get("provenance", {})
    provenance_rows = [
        [key, value]
        for key, value in sorted(provenance.items())
        if key != "excluded_activity_libraries"
        and isinstance(value, (str, int, float, bool))
    ]
    library_rows = [
        [
            library["id"],
            library["context"],
            library["assay"],
            library["layout"],
            library["qc_status"],
            library.get("estimated_fragment_length_bp"),
            library.get("source_project", ""),
            library.get("source_run_id", ""),
        ]
        for library in activity["libraries"]
    ] + [
        [
            library["id"],
            library["context"],
            library["assay"],
            library["layout"],
            library["qc_status"],
            library.get("estimated_fragment_length_bp"),
            "",
            library.get("reason", ""),
        ]
        for library in excluded
    ]

    source_inventory = [
        [
            record["kind"],
            str(record["path"]),
            record["path"].stat().st_size,
            record["observed_sha256"][:16],
        ]
        for record in records
    ]
    output_inventory = [
        [name, str(path.resolve()), path.stat().st_size, sha256_file(path)[:16]]
        for name, path in sorted(current_files.items())
    ]

    tss_profiles = [record for record in records if record["kind"] == "tss_profile"]
    profile_html = "".join(
        '<figure class="profile">'
        f'<img src="{_data_uri(record["path"])}" alt="{html.escape(record["path"].stem)} TSS profile"/>'
        f'<figcaption>{html.escape(record["path"].name.removesuffix(".profile.png"))}</figcaption>'
        "</figure>"
        for record in tss_profiles
    ) or '<p class="muted">No upstream TSS-profile plots were supplied.</p>'

    master_context_peaks = [
        (context, float(master_metrics.get("context_peak_counts", {}).get(context, 0)))
        for context in contexts
    ]
    mapped_chart = _bar_chart(
        "Mapped alignments by accepted/reviewed library",
        [
            (str(row["sample"]), float(row.get("mapped_alignments", 0)))
            for row in qc_rows
        ],
        color="#3976a8",
    )
    frip_chart = _bar_chart(
        "FRiP by accepted/reviewed library",
        [(str(row["sample"]), float(row.get("frip", 0))) for row in qc_rows],
        value_format="fraction",
        color="#3b8f6f",
    )
    factor_chart = _bar_chart(
        "Background-TMM normalization factors",
        [
            (f"{row['library_id']} ({row['assay']})", float(row["tmm_normalization_factor"]))
            for row in factors
        ],
        color="#df8f44",
    )
    nearest_active_counts = nearest_tss_metrics.get(
        "enhancer_with_active_nearest_promoter_by_context", {}
    )
    nearest_tss_section = f"""
<section class="stage"><h2>6. Default nearest-TSS enhancer candidates</h2>
{_cards([('Enhancers', nearest_tss_metrics.get('enhancer_count')), ('Enhancers with tied nearest TSSs', nearest_tss_metrics.get('nearest_tss_tie_element_count')), ('No same-chromosome TSS', nearest_tss_metrics.get('no_same_chromosome_promoter_count')), ('Assignment', 'nearest annotated TSS')])}
<p>This contact-independent branch has one row per proximal or distal enhancer. It retains every exact nearest-distance TSS tie in that row, reports gene and transcript identifiers, and adds per-context element activity plus whether any tied nearest promoter is active. Promoter activity requires a catalog-active context-member master-DHS summit within inclusive ±500 bp of the exact TSS. Enhancers are not filtered from this atlas-wide table by their own activity state.</p>
{_bar_chart('Enhancers whose nearest promoter is active by context', [(context, float(nearest_active_counts.get(context, 0))) for context in contexts], color='#3b8f6f')}
</section>
"""

    contact_section = ""
    if contact_metrics is not None:
        contact_rows = [
            [
                context,
                values["contact_strategy"],
                values["contact_assay"],
                values["contact_match"],
                values["contact_resolution_bp"],
                values["active_promoter_count"],
                values["element_promoter_edge_count"],
                values["element_gene_candidate_count"],
                values.get("active_contact_enhancer_gene_candidate_count", 0),
                values.get("nearest_active_promoter_gene_candidate_count", 0),
                values.get("active_distance_enhancer_gene_candidate_count", 0),
            ]
            for context, values in sorted(contact_metrics["contexts"].items())
        ]
        contact_section = f"""
<section class="stage"><h2>7. Optional context contact links and candidate genes</h2>
{_cards([('Observed contact contexts', contact_metrics.get('observed_context_count')), ('Distance-model contexts', contact_metrics.get('powerlaw_context_count')), ('Element–promoter edges', contact_metrics.get('element_promoter_edge_count')), ('Element–gene candidates', contact_metrics.get('element_gene_candidate_count')), ('Focused observed-contact candidates', contact_metrics.get('active_contact_enhancer_gene_candidate_count')), ('Nearest-active-promoter candidates', contact_metrics.get('nearest_active_promoter_gene_candidate_count')), ('Focused distance candidates', contact_metrics.get('active_distance_enhancer_gene_candidate_count'))])}
<p>This opt-in branch uses merged, ICE-balanced contact matrices for observed contexts. Configured contexts without a defensible map use the atlas-wide contact-decay model and remain explicitly labeled as distance-model evidence. Candidate scores combine contact weight with the promoter's context-resolved ATAC/H3K27ac activity. The per-context nearest-active-TSS projection inside this branch is supplementary; the wide default table above is independent of contact processing. Focused distance tables remain power-law rankings rather than measured contacts.</p>
{_table(['Context', 'Strategy', 'Assay', 'Match', 'Resolution bp', 'Active promoters', 'Element–promoter edges', 'Element–gene candidates', 'Focused observed', 'Nearest active TSS', 'Focused distance'], contact_rows, compact=True)}
</section>
"""

    css = """
    @page { size: A4 landscape; margin: 12mm 13mm 14mm; @bottom-right { content: "Page " counter(page) " of " counter(pages); color: #657184; font-size: 8pt; } }
    * { box-sizing: border-box; }
    body { color:#253047; font: 9.2pt/1.35 Arial, sans-serif; margin:0 auto; max-width:1120px; }
    h1 { color:#173f61; font-size:23pt; margin:0 0 4px; }
    h2 { color:#225a7f; border-bottom:2px solid #dbe6ef; padding-bottom:4px; margin:22px 0 10px; }
    h3 { color:#365b76; margin:16px 0 7px; }
    p { margin:5px 0 9px; }
    .subtitle { color:#607184; font-size:11pt; }
    .stage { page-break-before:always; }
    .cards { display:flex; flex-wrap:wrap; gap:8px; margin:10px 0 14px; }
    .card { background:#eef4f8; border-left:4px solid #3976a8; padding:8px 12px; min-width:145px; }
    .card-value { font-size:16pt; font-weight:700; color:#173f61; }
    .card-label { color:#607184; font-size:8pt; }
    table { border-collapse:collapse; width:100%; margin:7px 0 13px; font-size:8pt; }
    th { background:#e7eef4; color:#263c50; text-align:left; }
    th, td { border:1px solid #ccd7e1; padding:4px 5px; vertical-align:top; overflow-wrap:anywhere; }
    tr:nth-child(even) td { background:#f8fafc; }
    table.compact { font-size:7pt; }
    .warning { background:#fff4dc; border-left:4px solid #d68a13; padding:7px 10px; margin:5px 0; }
    .ok { background:#e8f5ee; border-left:4px solid #2e8b62; padding:7px 10px; }
    .muted { color:#738091; }
    .chart { width:100%; max-height:720px; margin:4px 0 12px; }
    .chart text { font:10px Arial,sans-serif; fill:#2b3544; }
    .chart .chart-title { font-size:13px; font-weight:bold; fill:#173f61; }
    .profiles { display:grid; grid-template-columns:repeat(4, 1fr); gap:7px; }
    .profile { margin:0; border:1px solid #d6dfe7; padding:3px; break-inside:avoid; }
    .profile img { width:100%; display:block; }
    .profile figcaption { text-align:center; color:#59697b; font-size:7pt; }
    .two-column { display:grid; grid-template-columns:1fr 1fr; gap:13px; }
    code { font-size:7.5pt; }
    """
    warning_html = (
        "".join(f'<div class="warning">{html.escape(value)}</div>' for value in warnings)
        if warnings
        else '<div class="ok">No rejected libraries or missing upstream QC rows were recorded.</div>'
    )
    report_html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{html.escape(config['project'])} integrated QC report</title><style>{css}</style></head>
<body>
<h1>Integrated atlas regulatory-catalog QC report</h1>
<p class="subtitle">Project <strong>{html.escape(str(config['project']))}</strong> · run <strong>{html.escape(str(config['run_id']))}</strong> · genome <strong>{html.escape(str(config['reference']['name']))}</strong></p>
{_cards([('Contexts', len(contexts)), ('Accepted libraries', len(activity['libraries'])), ('Excluded libraries', len(excluded)), ('Master DHSs', master_metrics.get('master_dhs_count')), ('Active calls', sum(sum(value.values()) for value in active_counts.values()))])}
<h2>Executive QC status</h2>{warning_html}
<p>The scientific endpoint uses exact variable-width master DHS ATAC signal, the maximum normalized density among three summit-centered H3K27ac windows, assay-specific 10-kb-background TMM, and a context-specific empirical background percentile. Activity requires a blacklist-clear context member above the 0.60 percentile of blacklist-clear DHSs that do not contribute to that context.</p>

<section class="stage"><h2>1. Inputs and provenance</h2>
<h3>Workflow provenance</h3>{_table(['Field', 'Value'], provenance_rows)}
<h3>Accepted and excluded libraries</h3>{_table(['Library', 'Context', 'Assay', 'Layout', 'QC status', 'SE fragment bp', 'Source project', 'Source run / reason'], library_rows, compact=True)}
<h3>Frozen upstream report inputs</h3>{_table(['Kind', 'Path', 'Bytes', 'SHA-256 prefix'], source_inventory, compact=True)}
</section>

<section class="stage"><h2>2. Read processing and library QC</h2>
{_cards([('QC libraries found', len(qc_rows)), ('Median mapped alignments', sorted([int(row.get('mapped_alignments', 0)) for row in qc_rows])[len(qc_rows)//2] if qc_rows else None), ('ATAC libraries', sum(row.get('assay') == 'atac' for row in qc_rows)), ('H3K27ac libraries', sum(row.get('assay') == 'chip_histone' for row in qc_rows))])}
<div class="two-column"><div>{mapped_chart}</div><div>{frip_chart}</div></div>
{_table(['Library', 'Context', 'Assay', 'Layout', 'Mapped alignments', 'FRiP', 'Peak count', 'Status'], qc_table_rows, compact=True)}
<h3>H3K27ac cross-correlation</h3>{_table(['Library', 'Reads', 'Estimated fragment shifts', 'NSC', 'RSC', 'Quality tag'], cross_correlation, compact=True)}
<h3>ATAC TSS profiles</h3><div class="profiles">{profile_html}</div>
</section>

<section class="stage"><h2>3. Replicate-supported master DHS registry</h2>
{_cards([('Source peaks', master_metrics.get('source_peak_count')), ('Master DHSs', master_metrics.get('master_dhs_count')), ('Multi-context DHSs', master_metrics.get('multi_context_master_dhs_count')), ('Median width bp', master_metrics.get('master_width_median')), ('Width range bp', f"{master_metrics.get('master_width_min')}–{master_metrics.get('master_width_max')}")])}
{_bar_chart('Replicate-supported peaks contributing by context', master_context_peaks, color='#7851a9')}
{_table(['Parameter', 'Value'], [[key, value] for key, value in sorted(master_metrics.items()) if key not in {'context_peak_counts', 'contexts'}])}
</section>

<section class="stage"><h2>4. Master-element quantification and background TMM</h2>
{_cards([('Master DHSs', activity_metrics.get('master_dhs_count')), ('Contexts', activity_metrics.get('context_count')), ('Libraries', activity_metrics.get('library_count')), ('Method', activity_metrics.get('normalization_method'))])}
{factor_chart}
{_table(['Library', 'Context', 'Assay', 'Usable units', 'TMM factor', 'Effective library size', 'Background features'], factor_rows, compact=True)}
</section>

<section class="stage"><h2>5. Regulatory-element catalog</h2>
{_cards([('Catalog rows', catalog_metrics.get('catalog_row_count')), ('Master DHSs', catalog_metrics.get('master_dhs_count')), ('Contexts', catalog_metrics.get('context_count')), ('Activity threshold', '> 0.60 background percentile')])}
{_stacked_chart('Active elements by TSS-distance class', contexts, active_counts)}
</section>

{nearest_tss_section}

{contact_section}

<section class="stage"><h2>{'8' if contact_metrics is not None else '7'}. Outputs and audit trail</h2>
<p>Every listed artifact is an explicit workflow input to this report. Sizes and checksum prefixes make the report auditable without embedding the large count tables themselves.</p>
{_table(['Output', 'Path', 'Bytes', 'SHA-256 prefix'], output_inventory, compact=True)}
</section>
</body></html>"""

    _atomic_text_if_changed(output_html, report_html)
    renderer = pdf_renderer or _default_pdf_renderer
    renderer_version = renderer(report_html, output_pdf)
    metrics = {
        "status": "ok",
        "schema_version": 2,
        "method": "integrated_regulatory_catalog_links_qc_report_v6",
        "project": config["project"],
        "run_id": config["run_id"],
        "contexts": contexts,
        "source_file_count": len(records),
        "current_artifact_count": len(current_files),
        "accepted_library_count": len(activity["libraries"]),
        "excluded_library_count": len(excluded),
        "upstream_qc_library_count": len(qc_rows),
        "nearest_tss_links": nearest_tss_metrics,
        "contact_links": contact_metrics,
        "warnings": warnings,
        "pdf_renderer": {"implementation": "WeasyPrint", "version": renderer_version},
        "inputs": {
            str(record["path"]): record["observed_sha256"] for record in records
        },
        "current_artifacts": {
            name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
            for name, path in sorted(current_files.items())
        },
        "outputs": {
            "html": {
                "path": str(output_html.resolve()),
                "sha256": sha256_file(output_html),
            },
            "pdf": {
                "path": str(output_pdf.resolve()),
                "sha256": sha256_file(output_pdf),
            },
        },
    }
    write_json_if_changed(output_metrics, metrics)
    return metrics
