import csv
from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.configuration import _reference_config
from short_read_processing.contact_metadata import (
    DM6_ATLAS_CONTACT_CONTEXTS,
    DM6_ATLAS_CONTEXT_IDS,
    default_dm6_contact_config,
    read_contact_source_manifest,
    verify_reported_checksum,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_MANIFEST = REPO_ROOT / "resources" / "atlas_contact_sources.tsv"


def test_canonical_contact_manifest_covers_seven_observed_contexts():
    rows = read_contact_source_manifest(SOURCE_MANIFEST)

    assert {row["context"] for row in rows} == {
        row["id"]
        for row in DM6_ATLAS_CONTACT_CONTEXTS
        if row["strategy"] == "observed"
    }
    assert len(rows) == 16
    assert {row["format"] for row in rows} == {"mcool", "cool.gz", "h5"}


def test_default_contact_config_models_exactly_two_contexts(tmp_path):
    reference = _reference_config("dm6", tmp_path / "references", tmp_path)

    config = default_dm6_contact_config(
        contexts=sorted(DM6_ATLAS_CONTEXT_IDS),
        reference=reference,
        manifest_path=SOURCE_MANIFEST,
        path_base=tmp_path,
    )

    assert config is not None
    assert [
        row["id"] for row in config["contexts"] if row["strategy"] == "powerlaw"
    ] == ["e13", "hid"]
    assert config["normalization"] == "merge_counts_then_ice_retry_v2"
    assert config["promoter_width_bp"] == 1000
    assert config["promoter_activity"] == "summit_within_500bp_catalog_activity_v3"
    assert config["maximum_distance_bp"] == 1_000_000
    assert config["candidate_observed_over_expected_threshold"] == 1.0
    assert len(config["source_manifest_sha256"]) == 64


def test_default_contact_config_is_not_applied_to_partial_context_sets(tmp_path):
    reference = _reference_config("dm6", tmp_path / "references", tmp_path)

    assert default_dm6_contact_config(
        contexts=["e5"],
        reference=reference,
        manifest_path=SOURCE_MANIFEST,
        path_base=tmp_path,
    ) is None


def test_contact_manifest_rejects_invalid_reported_checksum(tmp_path):
    rows = read_contact_source_manifest(SOURCE_MANIFEST)
    rows[0]["checksum"] = "sha256:not-a-digest"
    path = tmp_path / "contacts.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    with pytest.raises(AcquisitionError, match="invalid checksum"):
        read_contact_source_manifest(path)


def test_reported_checksum_is_checked_for_an_existing_source(tmp_path):
    source = tmp_path / "source.cool"
    source.write_bytes(b"contact matrix")

    with pytest.raises(AcquisitionError, match="checksum mismatch"):
        verify_reported_checksum(source, "sha256:" + "0" * 64)
