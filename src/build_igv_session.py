#!/usr/bin/env python3
"""Build portable IGV sessions for regulatory-atlas outputs."""

from __future__ import annotations

import argparse
import os
import xml.etree.ElementTree as ET
from pathlib import Path
import tempfile


def relative_path(path: Path, session_path: Path) -> str:
    return Path(os.path.relpath(path, session_path.parent)).as_posix()


def add_track(
    resources: ET.Element,
    panel: ET.Element,
    *,
    path: Path,
    output: Path,
    name: str,
    color: str,
    signal: bool,
) -> None:
    track_id = relative_path(path, output)
    ET.SubElement(resources, "Resource", path=track_id)
    ET.SubElement(
        panel,
        "Track",
        id=track_id,
        name=name,
        color=color,
        renderer="BAR_CHART" if signal else "BASIC_FEATURE",
        height="52" if signal else "24",
        expand="false",
        visible="true",
        windowFunction="mean" if signal else "count",
    )


def add_hidden_sequence_order_anchor(panel: ET.Element) -> None:
    """Reserve IGV's zero order so the master BED stays below genome tracks."""

    ET.SubElement(
        panel,
        "Track",
        id="Reference sequence",
        name="Reference sequence",
        clazz="org.broad.igv.track.SequenceTrack",
        visible="false",
    )


def _require(paths: list[Path]) -> None:
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing IGV inputs:\n" + "\n".join(map(str, missing)))


def _write_xml_if_changed(session: ET.Element, output: Path) -> None:
    ET.indent(session, space="  ")
    content = ET.tostring(session, encoding="utf-8", xml_declaration=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file() and output.read_bytes() == content:
        return
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{output.name}.",
            dir=output.parent,
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(content)
        os.replace(temporary_name, output)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def build_session(
    atac_root: Path,
    output: Path,
    genome: str,
    locus: str,
    *,
    master_bed: Path | None = None,
    catalog_bed_root: Path,
    catalog_track_root: Path | None = None,
) -> tuple[int, int]:
    """Build one master-first session with five ordered tracks per context."""

    conditions_root = atac_root / "conditions"
    conditions = sorted(path for path in conditions_root.iterdir() if path.is_dir())
    if not conditions:
        raise ValueError(f"No ATAC condition outputs found under {conditions_root}")
    catalog_track_root = catalog_track_root or catalog_bed_root.parent / "tracks"
    context_beds = {
        path.name.removesuffix(".dhs.bed"): path
        for path in catalog_bed_root.glob("*.dhs.bed")
    }
    condition_names = {path.name for path in conditions}
    if set(context_beds) != condition_names:
        missing_catalog = sorted(condition_names - set(context_beds))
        missing_atac = sorted(set(context_beds) - condition_names)
        details = []
        if missing_catalog:
            details.append("missing catalog contexts: " + ", ".join(missing_catalog))
        if missing_atac:
            details.append("missing ATAC conditions: " + ", ".join(missing_atac))
        raise ValueError("ATAC and catalog contexts differ; " + "; ".join(details))

    output.parent.mkdir(parents=True, exist_ok=True)
    session = ET.Element(
        "Session",
        genome=genome,
        hasGeneTrack="true",
        hasSequenceTrack="false",
        locus=locus,
        version="3",
    )
    resources = ET.SubElement(session, "Resources")
    panel = ET.SubElement(session, "Panel", name="Context regulatory atlas")
    add_hidden_sequence_order_anchor(panel)
    if master_bed is None:
        master_bed = atac_root / "master" / "master_dhs.bed"
    _require([master_bed])
    add_track(
        resources,
        panel,
        path=master_bed,
        output=output,
        name="Master DHS registry",
        color="106,27,154",
        signal=False,
    )
    track_count = 1
    for condition_root in conditions:
        condition = condition_root.name
        label = condition.upper()
        tracks = condition_root / "tracks"
        context_dhs = context_beds[condition]
        inputs = [
            catalog_track_root
            / f"{condition}.atac.mean.background_tmm.bw",
            tracks / f"{condition}.qpois.bw",
            catalog_track_root
            / f"{condition}.h3k27ac.mean.background_tmm.bw",
            context_dhs,
            catalog_bed_root / f"{condition}.elements.bed",
        ]
        _require(inputs)
        specifications = (
            (
                inputs[0],
                f"{label} | mean ATAC 150-bp Tn5 pileup (background-TMM)",
                "31,120,180",
                True,
            ),
            (inputs[1], f"{label} | pooled ATAC qpois signal", "117,112,179", True),
            (
                inputs[2],
                f"{label} | mean H3K27ac signal (background-TMM)",
                "221,126,32",
                True,
            ),
            (inputs[3], f"{label} | context DHSs", "0,145,130", False),
            (inputs[4], f"{label} | candidate cCREs (activity-scored)", "202,61,52", False),
        )
        for path, name, color, signal in specifications:
            add_track(
                resources,
                panel,
                path=path,
                output=output,
                name=name,
                color=color,
                signal=signal,
            )
            track_count += 1

    _write_xml_if_changed(session, output)
    return len(conditions), track_count


def build_catalog_session(
    *,
    context: str,
    genome: str,
    atac_bigwig: Path,
    h3k27ac_bigwig: Path,
    context_dhs_bed: Path,
    master_dhs_bed: Path,
    elements_bed: Path,
    output: Path,
    locus: str = "All",
) -> int:
    """Build a portable five-track IGV session for one atlas context."""

    paths = [
        atac_bigwig,
        h3k27ac_bigwig,
        context_dhs_bed,
        master_dhs_bed,
        elements_bed,
    ]
    _require(paths)
    session = ET.Element(
        "Session",
        genome=genome,
        hasGeneTrack="true",
        hasSequenceTrack="false",
        locus=locus,
        version="3",
    )
    resources = ET.SubElement(session, "Resources")
    panel = ET.SubElement(session, "Panel", name=f"{context} regulatory atlas")
    add_hidden_sequence_order_anchor(panel)
    label = context.upper()
    specifications = [
        (
            master_dhs_bed,
            "Master DHS registry",
            "106,27,154",
            False,
        ),
        (
            atac_bigwig,
            f"{label} | mean ATAC 150-bp Tn5 pileup (background-TMM)",
            "31,120,180",
            True,
        ),
        (
            h3k27ac_bigwig,
            f"{label} | mean H3K27ac fragment coverage (background-TMM)",
            "221,126,32",
            True,
        ),
        (
            context_dhs_bed,
            f"{label} | context DHSs (master coordinates)",
            "0,145,130",
            False,
        ),
        (
            elements_bed,
            f"{label} | candidate cCREs (activity-scored)",
            "202,61,52",
            False,
        ),
    ]
    for path, name, color, signal in specifications:
        add_track(
            resources,
            panel,
            path=path,
            output=output,
            name=name,
            color=color,
            signal=signal,
        )
    _write_xml_if_changed(session, output)
    return len(specifications)


def build_all_contexts_catalog_session(
    *,
    contexts: list[str],
    genome: str,
    atac_bigwigs: dict[str, Path],
    h3k27ac_bigwigs: dict[str, Path],
    context_dhs_beds: dict[str, Path],
    master_dhs_bed: Path,
    element_beds: dict[str, Path],
    output: Path,
    locus: str = "All",
) -> int:
    """Build one portable IGV session containing every catalog context."""

    if not contexts or len(contexts) != len(set(contexts)):
        raise ValueError("IGV catalog contexts must be non-empty and unique")
    mappings = {
        "ATAC BigWigs": atac_bigwigs,
        "H3K27ac BigWigs": h3k27ac_bigwigs,
        "context DHS BEDs": context_dhs_beds,
        "activity-annotated element BEDs": element_beds,
    }
    expected = set(contexts)
    for label, paths in mappings.items():
        if set(paths) != expected:
            raise ValueError(f"{label} must cover exactly the IGV catalog contexts")
    _require(
        [
            master_dhs_bed,
            *(path for paths in mappings.values() for path in paths.values()),
        ]
    )

    session = ET.Element(
        "Session",
        genome=genome,
        hasGeneTrack="true",
        hasSequenceTrack="false",
        locus=locus,
        version="3",
    )
    resources = ET.SubElement(session, "Resources")
    panel = ET.SubElement(session, "Panel", name="All-context regulatory atlas")
    add_hidden_sequence_order_anchor(panel)
    add_track(
        resources,
        panel,
        path=master_dhs_bed,
        output=output,
        name="Master DHS registry",
        color="106,27,154",
        signal=False,
    )
    track_count = 1
    for context in contexts:
        label = context.upper()
        for path, name, color, signal in (
            (
                atac_bigwigs[context],
                f"{label} | mean ATAC 150-bp Tn5 pileup (background-TMM)",
                "31,120,180",
                True,
            ),
            (
                h3k27ac_bigwigs[context],
                f"{label} | mean H3K27ac fragment coverage (background-TMM)",
                "221,126,32",
                True,
            ),
            (
                context_dhs_beds[context],
                f"{label} | context DHSs (master coordinates)",
                "0,145,130",
                False,
            ),
            (
                element_beds[context],
                f"{label} | candidate cCREs (activity-scored)",
                "202,61,52",
                False,
            ),
        ):
            add_track(
                resources,
                panel,
                path=path,
                output=output,
                name=name,
                color=color,
                signal=signal,
            )
            track_count += 1
    _write_xml_if_changed(session, output)
    return track_count



def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("atac_root", type=Path, help="Run's results/.../atac directory")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genome", default="dm6")
    parser.add_argument("--locus", default="All")
    parser.add_argument(
        "--master-bed",
        type=Path,
        help="Master DHS BED; defaults to ATAC_ROOT/master/master_dhs.bed when present",
    )
    parser.add_argument(
        "--catalog-bed-root",
        type=Path,
        required=True,
        help="activity/catalog/bed directory with all context and activity-annotated BEDs",
    )
    parser.add_argument(
        "--catalog-track-root",
        type=Path,
        help="Directory with mean H3K27ac BigWigs; defaults to BED_ROOT/../tracks",
    )
    args = parser.parse_args()
    condition_n, track_n = build_session(
        args.atac_root.resolve(),
        args.output.resolve(),
        args.genome,
        args.locus,
        master_bed=args.master_bed.resolve() if args.master_bed else None,
        catalog_bed_root=args.catalog_bed_root.resolve(),
        catalog_track_root=(
            args.catalog_track_root.resolve() if args.catalog_track_root else None
        ),
    )
    print(
        f"Wrote {condition_n} contexts and {track_n} tracks to "
        f"{args.output.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
