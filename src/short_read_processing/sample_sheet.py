"""Parse and validate the accession-first CSV/TSV sample sheet."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import yaml

from .accessions import AcquisitionError, normalize_accession


DEFAULT_SCHEMA = Path(__file__).resolve().parents[2] / "schemas" / "sample-sheet.schema.yaml"
TECHNICAL_RUN_FIELDS = {"accession", "notes"}
ASSAY_ALIASES = {"h3k27ac": "chip_histone"}
MACS3_FIELDS = {
    "macs3_format",
    "macs3_qvalue",
    "macs3_broad",
    "macs3_broad_cutoff",
    "macs3_nomodel",
    "macs3_shift",
    "macs3_extsize",
}


def delimiter_for(path: Path, sample: str) -> str:
    """Choose comma or tab deterministically, using the suffix when available."""

    suffix = path.suffix.lower()
    if suffix == ".csv":
        return ","
    if suffix in {".tsv", ".tab"}:
        return "\t"
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error as exc:
        raise AcquisitionError(f"Could not determine whether {path} is CSV or TSV") from exc


def read_delimited_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a UTF-8 CSV/TSV and return stripped headers and rows."""

    text = path.read_text(encoding="utf-8-sig")
    delimiter = delimiter_for(path, text[:8192])
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    fieldnames = [field.strip() for field in (reader.fieldnames or []) if field is not None]
    if not fieldnames:
        raise AcquisitionError(f"Sample sheet {path} has no header")
    if len(fieldnames) != len(set(fieldnames)):
        raise AcquisitionError(f"Sample sheet {path} has duplicate column names")
    rows: list[dict[str, str]] = []
    for row in reader:
        normalized = {
            str(key).strip(): (value or "").strip()
            for key, value in row.items()
            if key is not None
        }
        if any(normalized.values()):
            rows.append(normalized)
    return fieldnames, rows


def load_schema(path: Path = DEFAULT_SCHEMA) -> dict[str, Any]:
    schema = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(schema, dict) or not isinstance(schema.get("fields"), dict):
        raise AcquisitionError(f"Invalid sample-sheet schema: {path}")
    return schema


def _parse_boolean(value: str, *, field: str, line: int) -> bool:
    normalized = value.lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise AcquisitionError(f"Line {line}: {field} must be true or false")


def _parse_field(name: str, value: str, spec: dict[str, Any], *, line: int) -> Any:
    if not value:
        if spec.get("allow_blank"):
            return None
        if "default" in spec:
            return spec["default"]
        raise AcquisitionError(f"Line {line}: required value {name!r} is blank")

    field_type = spec.get("type", "string")
    try:
        if field_type == "accession":
            parsed: Any = normalize_accession(value)
        elif field_type == "integer":
            if not re.fullmatch(r"[+-]?\d+", value):
                raise ValueError
            parsed = int(value)
        elif field_type == "number":
            parsed = float(value)
        elif field_type == "boolean":
            parsed = _parse_boolean(value, field=name, line=line)
        else:
            parsed = value
    except ValueError as exc:
        raise AcquisitionError(f"Line {line}: invalid {field_type} value for {name}: {value!r}") from exc

    if field_type == "enum" and parsed not in spec.get("values", []):
        choices = ", ".join(str(item) for item in spec.get("values", []))
        raise AcquisitionError(f"Line {line}: {name} must be one of {choices}")
    if field_type in {"string", "accession"} and spec.get("pattern"):
        if not re.fullmatch(str(spec["pattern"]), str(parsed)):
            raise AcquisitionError(f"Line {line}: invalid {name}: {parsed!r}")
    if isinstance(parsed, (int, float)) and not isinstance(parsed, bool):
        if "minimum" in spec and parsed < spec["minimum"]:
            raise AcquisitionError(f"Line {line}: {name} must be at least {spec['minimum']}")
        if "exclusive_minimum" in spec and parsed <= spec["exclusive_minimum"]:
            raise AcquisitionError(f"Line {line}: {name} must be greater than {spec['exclusive_minimum']}")
        if "maximum" in spec and parsed > spec["maximum"]:
            raise AcquisitionError(f"Line {line}: {name} must be at most {spec['maximum']}")
    return parsed


def _resolve_defaults(row: dict[str, Any], schema: dict[str, Any], *, line: int) -> None:
    assay = str(row["assay"])
    fields = schema["fields"]
    supplied = row["__supplied__"]
    for name, spec in fields.items():
        if (
            row.get(name) is None
            and "defaults_by_assay" in spec
            and name not in MACS3_FIELDS
        ):
            row[name] = spec["defaults_by_assay"].get(assay)
        if (
            row.get(name) is None
            and "default" in spec
            and name not in MACS3_FIELDS
        ):
            row[name] = spec["default"]

    for name in MACS3_FIELDS:
        if row.get(name) is None and "defaults_by_assay" in fields[name]:
            row[name] = fields[name]["defaults_by_assay"].get(assay)
        if row.get(name) is None and "default" in fields[name]:
            row[name] = fields[name]["default"]
    if assay == "chip_histone" and not row["macs3_broad"]:
        raise AcquisitionError(
            f"Line {line}: H3K27ac peaks are called with macs3_broad=true"
        )
    if not row["macs3_broad"]:
        if "macs3_broad_cutoff" in supplied:
            raise AcquisitionError(
                f"Line {line}: macs3_broad_cutoff requires macs3_broad=true"
            )
        row["macs3_broad_cutoff"] = None
    if row["macs3_shift"] is not None or row["macs3_extsize"] is not None:
        if assay != "atac" and row["macs3_format"] != "BAM":
            raise AcquisitionError(
                f"Line {line}: H3K27ac macs3_shift/extsize requires macs3_format=BAM"
            )
        if row["macs3_extsize"] is None:
            raise AcquisitionError(f"Line {line}: macs3_shift requires macs3_extsize")
        if "macs3_nomodel" in supplied and not row["macs3_nomodel"]:
            raise AcquisitionError(
                f"Line {line}: macs3_shift/extsize conflicts with macs3_nomodel=false"
            )
        row["macs3_nomodel"] = True

    if row["adapter_preset"] == "custom" and not row.get("adapter_fasta"):
        raise AcquisitionError(f"Line {line}: adapter_preset=custom requires adapter_fasta")
    if row.get("adapter_fasta") and row["adapter_preset"] != "custom":
        raise AcquisitionError(f"Line {line}: adapter_fasta requires adapter_preset=custom")


def _validate_relationships(rows: list[dict[str, Any]], canonical_fields: set[str]) -> None:
    by_library: dict[str, list[dict[str, Any]]] = {}
    accession_lines: dict[str, int] = {}
    for row in rows:
        accession = str(row["accession"])
        if accession in accession_lines:
            raise AcquisitionError(
                f"Line {row['__line__']}: duplicate accession {accession}; "
                f"first seen on line {accession_lines[accession]}"
            )
        accession_lines[accession] = int(row["__line__"])
        by_library.setdefault(str(row["library_id"]), []).append(row)

    for library_id, library_rows in by_library.items():
        first = library_rows[0]
        for row in library_rows[1:]:
            mismatches = [
                field
                for field in canonical_fields
                if field not in TECHNICAL_RUN_FIELDS
                and row.get(field) != first.get(field)
            ]
            if mismatches:
                raise AcquisitionError(
                    f"Rows for library_id {library_id!r} disagree on: "
                    + ", ".join(mismatches)
                )


def read_sample_sheet(
    path: Path,
    *,
    schema_path: Path = DEFAULT_SCHEMA,
) -> list[dict[str, Any]]:
    """Return typed, default-resolved rows from the canonical sample sheet."""

    schema = load_schema(schema_path)
    fieldnames, raw_rows = read_delimited_rows(path)
    missing = [name for name in schema.get("required_columns", []) if name not in fieldnames]
    if missing:
        raise AcquisitionError(
            f"Sample sheet {path} is missing required columns: " + ", ".join(missing)
        )
    if not raw_rows:
        raise AcquisitionError(f"Sample sheet {path} has no data rows")
    unsupported = [name for name in ("role", "control_library") if name in fieldnames]
    if unsupported:
        raise AcquisitionError(
            f"Sample sheet {path} uses unsupported columns: "
            + ", ".join(unsupported)
            + "; H3K27ac libraries are processed IP-only"
        )

    fields = schema["fields"]
    parsed_rows: list[dict[str, Any]] = []
    for line, raw in enumerate(raw_rows, start=2):
        parsed = {
            name: _parse_field(name, raw.get(name, ""), spec, line=line)
            for name, spec in fields.items()
            if name in fieldnames
        }
        for name in fields:
            if name not in parsed:
                parsed[name] = None
        for name, value in raw.items():
            if name not in fields:
                parsed[name] = value
        parsed["__line__"] = line
        parsed["__supplied__"] = {name for name, value in raw.items() if value}
        parsed["assay"] = ASSAY_ALIASES.get(str(parsed["assay"]), parsed["assay"])
        parsed["role"] = "treatment"
        _resolve_defaults(parsed, schema, line=line)
        parsed_rows.append(parsed)

    _validate_relationships(parsed_rows, set(fields))
    return parsed_rows


def sample_sheet_accessions(path: Path, *, schema_path: Path = DEFAULT_SCHEMA) -> list[str]:
    rows = read_sample_sheet(path, schema_path=schema_path)
    return list(dict.fromkeys(str(row["accession"]) for row in rows))
