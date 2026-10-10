import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from build_igv_session import (
    build_all_contexts_catalog_session,
    build_catalog_session,
    build_session,
)


def touch_all(paths):
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()


def user_tracks(root):
    tracks = root.findall("./Panel/Track")
    assert tracks[0].attrib == {
        "id": "Reference sequence",
        "name": "Reference sequence",
        "clazz": "org.broad.igv.track.SequenceTrack",
        "visible": "false",
    }
    return tracks[1:]


def test_session_has_master_then_five_ordered_tracks_per_context(tmp_path: Path):
    atac = tmp_path / "atac"
    catalog = tmp_path / "catalog"
    contexts = ["e13", "e5"]
    paths = [atac / "master/master_dhs.bed"]
    for context in contexts:
        paths.extend(
            [
                catalog
                / f"tracks/{context}.atac.mean.background_tmm.bw",
                atac / f"conditions/{context}/tracks/{context}.qpois.bw",
                catalog
                / f"tracks/{context}.h3k27ac.mean.background_tmm.bw",
                catalog / f"bed/{context}.dhs.bed",
                catalog / f"bed/{context}.elements.bed",
            ]
        )
    touch_all(paths)
    output = tmp_path / "session.xml"

    counts = build_session(
        atac,
        output,
        "dm6",
        "All",
        catalog_bed_root=catalog / "bed",
    )

    root = ET.parse(output).getroot()
    assert root.attrib["hasGeneTrack"] == "true"
    assert root.attrib["hasSequenceTrack"] == "false"
    resources = root.findall("./Resources/Resource")
    tracks = user_tracks(root)
    assert counts == (2, 11)
    assert len(resources) == len(tracks) == 11
    assert resources[0].attrib["path"] == "atac/master/master_dhs.bed"
    assert all((output.parent / item.attrib["path"]).is_file() for item in resources)
    assert [track.attrib["name"] for track in tracks] == [
        "Master DHS registry",
        "E13 | mean ATAC 150-bp Tn5 pileup (background-TMM)",
        "E13 | pooled ATAC qpois signal",
        "E13 | mean H3K27ac signal (background-TMM)",
        "E13 | context DHSs",
        "E13 | candidate cCREs (activity-scored)",
        "E5 | mean ATAC 150-bp Tn5 pileup (background-TMM)",
        "E5 | pooled ATAC qpois signal",
        "E5 | mean H3K27ac signal (background-TMM)",
        "E5 | context DHSs",
        "E5 | candidate cCREs (activity-scored)",
    ]


def test_session_requires_matching_atac_and_catalog_contexts(tmp_path: Path):
    atac = tmp_path / "atac"
    condition = atac / "conditions" / "e5"
    touch_all(
        [
            atac / "master/master_dhs.bed",
            condition / "tracks/e5.MACS3-pileup.unscaled.bw",
            condition / "tracks/e5.qpois.bw",
            tmp_path / "catalog/bed/e13.dhs.bed",
            tmp_path / "catalog/bed/e13.elements.bed",
        ]
    )
    output = tmp_path / "session.xml"

    with pytest.raises(ValueError, match="ATAC and catalog contexts differ"):
        build_session(
            atac,
            output,
            "dm6",
            "All",
            catalog_bed_root=tmp_path / "catalog/bed",
        )


def test_catalog_session_has_relative_mean_signal_and_element_tracks(tmp_path: Path):
    inputs = {
        "atac_bigwig": tmp_path / "catalog/tracks/ctx.atac.bw",
        "h3k27ac_bigwig": tmp_path / "catalog/tracks/ctx.h3k27ac.bw",
        "context_dhs_bed": tmp_path / "catalog/bed/ctx.dhs.bed",
        "master_dhs_bed": tmp_path / "catalog/bed/master_dhs.bed",
        "elements_bed": tmp_path / "catalog/bed/ctx.elements.bed",
    }
    touch_all(list(inputs.values()))
    output = tmp_path / "catalog/igv/ctx.xml"

    count = build_catalog_session(
        context="ctx",
        genome="dm6",
        output=output,
        **inputs,
    )

    root = ET.parse(output).getroot()
    assert root.attrib["hasGeneTrack"] == "true"
    assert root.attrib["hasSequenceTrack"] == "false"
    resources = root.findall("./Resources/Resource")
    tracks = user_tracks(root)
    assert count == len(resources) == len(tracks) == 5
    assert resources[0].attrib["path"] == "../bed/master_dhs.bed"
    assert all((output.parent / item.attrib["path"]).is_file() for item in resources)
    assert [track.attrib["name"] for track in tracks] == [
        "Master DHS registry",
        "CTX | mean ATAC 150-bp Tn5 pileup (background-TMM)",
        "CTX | mean H3K27ac fragment coverage (background-TMM)",
        "CTX | context DHSs (master coordinates)",
        "CTX | candidate cCREs (activity-scored)",
    ]


def test_all_contexts_catalog_session_contains_every_context_once(tmp_path: Path):
    contexts = ["ctx_a", "ctx_b"]
    master = tmp_path / "catalog/bed/master_dhs.bed"
    atac = {
        context: tmp_path / f"catalog/tracks/{context}.atac.bw"
        for context in contexts
    }
    h3k27ac = {
        context: tmp_path / f"catalog/tracks/{context}.h3k27ac.bw"
        for context in contexts
    }
    context_dhs = {
        context: tmp_path / f"catalog/bed/{context}.dhs.bed"
        for context in contexts
    }
    elements = {
        context: tmp_path / f"catalog/bed/{context}.elements.bed"
        for context in contexts
    }
    touch_all([master, *atac.values(), *h3k27ac.values(), *context_dhs.values(), *elements.values()])
    output = tmp_path / "catalog/all-contexts.igv.xml"

    count = build_all_contexts_catalog_session(
        contexts=contexts,
        genome="dm6",
        atac_bigwigs=atac,
        h3k27ac_bigwigs=h3k27ac,
        context_dhs_beds=context_dhs,
        master_dhs_bed=master,
        element_beds=elements,
        output=output,
    )

    root = ET.parse(output).getroot()
    assert root.attrib["hasGeneTrack"] == "true"
    assert root.attrib["hasSequenceTrack"] == "false"
    resources = root.findall("./Resources/Resource")
    tracks = user_tracks(root)
    assert count == len(resources) == len(tracks) == 9
    assert resources[0].attrib["path"] == "bed/master_dhs.bed"
    assert all((output.parent / item.attrib["path"]).is_file() for item in resources)
    assert [track.attrib["name"] for track in tracks] == [
        "Master DHS registry",
        "CTX_A | mean ATAC 150-bp Tn5 pileup (background-TMM)",
        "CTX_A | mean H3K27ac fragment coverage (background-TMM)",
        "CTX_A | context DHSs (master coordinates)",
        "CTX_A | candidate cCREs (activity-scored)",
        "CTX_B | mean ATAC 150-bp Tn5 pileup (background-TMM)",
        "CTX_B | mean H3K27ac fragment coverage (background-TMM)",
        "CTX_B | context DHSs (master coordinates)",
        "CTX_B | candidate cCREs (activity-scored)",
    ]


def test_all_contexts_catalog_session_requires_complete_context_mappings(tmp_path: Path):
    with pytest.raises(ValueError, match="ATAC BigWigs"):
        build_all_contexts_catalog_session(
            contexts=["ctx_a", "ctx_b"],
            genome="dm6",
            atac_bigwigs={"ctx_a": tmp_path / "ctx_a.atac.bw"},
            h3k27ac_bigwigs={
                context: tmp_path / f"{context}.h3k27ac.bw"
                for context in ("ctx_a", "ctx_b")
            },
            context_dhs_beds={
                context: tmp_path / f"{context}.dhs.bed"
                for context in ("ctx_a", "ctx_b")
            },
            master_dhs_bed=tmp_path / "master.bed",
            element_beds={
                context: tmp_path / f"{context}.elements.bed"
                for context in ("ctx_a", "ctx_b")
            },
            output=tmp_path / "all-contexts.igv.xml",
        )
