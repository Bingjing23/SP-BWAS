#!/usr/bin/env python3

"""
Summarize disease counts in usable rsfMRI participant cohorts.

Purpose:
    Intersect processed/usable rsfMRI EID lists with a phenotype table and
    report case/control/missing counts for disease groups used in a separate
    neuroimaging project. This script supports either one EID list or a cohort
    manifest with primary and supplementary cohorts.

Inputs:
    --phenotype-table: wide table with an EID column and disease phenotype
        columns, where cases are coded by values listed in the disease
        definition table.
    --usable-list: optional single-cohort file with one EID column.
    --cohort-manifest: optional multi-cohort table with cohort_id,
        cohort_label, usable_list, analysis_role, and notes.
    --definitions: disease grouping table; defaults to
        manifests/rsfmri_disease_definitions.tsv.

Outputs:
    outputs/rsfmri_disease_profile/rsfmri_disease_counts.tsv
    outputs/rsfmri_disease_profile/rsfmri_group_counts.tsv
    outputs/rsfmri_disease_profile/rsfmri_disease_missing_definitions.tsv
    outputs/rsfmri_disease_profile/subjects_with_disease_in_gwas_core.tsv
    outputs/rsfmri_disease_profile/subjects_with_disease_in_embedding_cohort.tsv
    outputs/rsfmri_disease_profile/disease_positive_embedding_not_gwas_core.tsv
    outputs/rsfmri_disease_profile/rsfmri_disease_summary_report.md

How to run:
    python3 scripts/08_summarize_rsfmri_disease_counts.py \
        --cohort-manifest manifests/rsfmri_cohorts.tsv \
        --phenotype-table /path/to/phenotype_table.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


EID_COLUMNS = (
    "eid",
    "EID",
    "f.eid",
    "f.eid.0.0",
    "IID",
    "id",
    "ID",
    "participant_id",
    "subject_id",
    "V1",
)


@dataclass(frozen=True)
class Cohort:
    cohort_id: str
    cohort_label: str
    usable_list: Path
    analysis_role: str
    notes: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize disease counts in usable rsfMRI participant cohorts."
    )
    parser.add_argument("--project-root", default=".", help="Project root directory.")
    parser.add_argument("--usable-list", default="", help="Single usable EID list.")
    parser.add_argument(
        "--cohort-manifest",
        default="",
        help="TSV manifest describing one or more usable EID cohorts.",
    )
    parser.add_argument(
        "--phenotype-table",
        required=True,
        help="Wide phenotype table with EID and disease phenotype columns.",
    )
    parser.add_argument(
        "--definitions",
        default="manifests/rsfmri_disease_definitions.tsv",
        help="Disease definition TSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs/rsfmri_disease_profile",
        help="Output directory.",
    )
    args = parser.parse_args()

    if bool(args.usable_list) == bool(args.cohort_manifest):
        parser.error("Provide exactly one of --usable-list or --cohort-manifest.")

    return args


def resolve_path(path: str | Path, root: Path) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return root / path


def detect_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        first_line = handle.readline()
    if "," in first_line and "\t" not in first_line:
        return ","
    return "\t"


def read_delimited(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Input file is empty: {path}")

    delimiter = detect_delimiter(path)
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Input file has no header: {path}")
        return list(reader)


def write_tsv(path: Path, rows: Sequence[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def normalize_eid(value: object) -> str:
    """Return bare participant EID from bare EID or case IDs like 1013356_20227_2_0."""
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN"}:
        return ""
    match = re.match(r"^(\d+)", text)
    if match:
        return match.group(1)
    return text


def read_usable_ids(path: Path) -> List[str]:
    if not path.exists():
        raise FileNotFoundError(f"Usable-list file does not exist: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Usable-list file is empty: {path}")

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        lines = [line.strip() for line in handle if line.strip()]

    if not lines:
        return []

    first_tokens = re.split(r"[\t, ]+", lines[0])
    first_token_eid = normalize_eid(first_tokens[0]) if first_tokens else ""
    has_header = (
        len(first_tokens) > 1 and any(token.lower() in {col.lower() for col in EID_COLUMNS} for token in first_tokens)
    ) or (
        len(first_tokens) == 1 and not first_token_eid.isdigit()
    )
    eid_col_index = 0

    if has_header:
        lowered = [token.lower() for token in first_tokens]
        for candidate in EID_COLUMNS:
            if candidate.lower() in lowered:
                eid_col_index = lowered.index(candidate.lower())
                break
        data_lines = lines[1:]
    else:
        data_lines = lines

    ids: List[str] = []
    seen = set()
    for line in data_lines:
        tokens = re.split(r"[\t, ]+", line)
        if eid_col_index >= len(tokens):
            continue
        eid = normalize_eid(tokens[eid_col_index])
        if eid and eid not in seen:
            seen.add(eid)
            ids.append(eid)

    return ids


def find_eid_column(fieldnames: Iterable[str], source: Path) -> str:
    fields = list(fieldnames)
    for candidate in EID_COLUMNS:
        if candidate in fields:
            return candidate

    lowered = {field.lower(): field for field in fields}
    for candidate in EID_COLUMNS:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]

    raise ValueError(f"Could not identify an EID column in {source}")


def load_phenotypes(path: Path) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    rows = read_delimited(path)
    if not rows:
        return {}, []

    fieldnames = list(rows[0].keys())
    eid_col = find_eid_column(fieldnames, path)
    phenotypes: Dict[str, Dict[str, str]] = {}
    duplicate_count = 0
    for row in rows:
        eid = normalize_eid(row.get(eid_col, ""))
        if not eid:
            continue
        if eid in phenotypes:
            duplicate_count += 1
            continue
        phenotypes[eid] = row

    if duplicate_count:
        print(f"Warning: ignored {duplicate_count} duplicate phenotype rows by EID.")
    return phenotypes, fieldnames


def split_tokens(value: str) -> List[str]:
    return [token.strip() for token in str(value).split(";") if token.strip()]


def load_definitions(path: Path) -> List[Dict[str, str]]:
    rows = read_delimited(path)
    required = {
        "disease_group",
        "disease_id",
        "disease_label",
        "candidate_columns",
        "positive_values",
        "notes",
    }
    missing = required.difference(rows[0].keys() if rows else set())
    if missing:
        raise ValueError(f"Definition table is missing columns: {', '.join(sorted(missing))}")
    return rows


def load_cohorts(args: argparse.Namespace, root: Path) -> List[Cohort]:
    if args.cohort_manifest:
        manifest_path = resolve_path(args.cohort_manifest, root)
        rows = read_delimited(manifest_path)
        required = {"cohort_id", "cohort_label", "usable_list", "analysis_role", "notes"}
        missing = required.difference(rows[0].keys() if rows else set())
        if missing:
            raise ValueError(f"Cohort manifest is missing columns: {', '.join(sorted(missing))}")
        return [
            Cohort(
                cohort_id=row["cohort_id"],
                cohort_label=row["cohort_label"],
                usable_list=resolve_path(row["usable_list"], root),
                analysis_role=row["analysis_role"],
                notes=row["notes"],
            )
            for row in rows
        ]

    return [
        Cohort(
            cohort_id="single_cohort",
            cohort_label="Single usable rsfMRI cohort",
            usable_list=resolve_path(args.usable_list, root),
            analysis_role="single",
            notes="",
        )
    ]


def find_phenotype_column(fieldnames: Sequence[str], candidates: Sequence[str]) -> str:
    field_set = set(fieldnames)
    for candidate in candidates:
        if candidate in field_set:
            return candidate

    lower_to_field = {field.lower(): field for field in fieldnames}
    for candidate in candidates:
        field = lower_to_field.get(candidate.lower())
        if field:
            return field

    return ""


def case_status(value: object, positive_values: Sequence[str]) -> str:
    text = str(value).strip()
    if text == "" or text.upper() in {"NA", "NAN", "NONE"}:
        return "missing"
    if text in set(positive_values):
        return "case"
    return "control"


def summarize_cohort(
    cohort: Cohort,
    phenotypes: Dict[str, Dict[str, str]],
    phenotype_columns: Sequence[str],
    definitions: Sequence[Dict[str, str]],
    out_dir: Path,
) -> Dict[str, object]:
    usable_ids = read_usable_ids(cohort.usable_list)
    usable_set = set(usable_ids)

    disease_counts: List[Dict[str, object]] = []
    subject_cases: List[Dict[str, object]] = []
    group_subject_status: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))

    for definition in definitions:
        candidates = split_tokens(definition["candidate_columns"])
        positive_values = split_tokens(definition["positive_values"])
        phenotype_column = find_phenotype_column(phenotype_columns, candidates)
        disease_group = definition["disease_group"]
        disease_id = definition["disease_id"]

        if not phenotype_column:
            disease_counts.append(
                {
                    "cohort_id": cohort.cohort_id,
                    "cohort_label": cohort.cohort_label,
                    "analysis_role": cohort.analysis_role,
                    "disease_group": disease_group,
                    "disease_id": disease_id,
                    "disease_label": definition["disease_label"],
                    "phenotype_column": "",
                    "usable_total_n": len(usable_ids),
                    "case_n": "",
                    "control_n": "",
                    "missing_n": "",
                    "nonmissing_n": "",
                    "prevalence": "",
                    "ready": "false",
                    "issue": "no matching phenotype column",
                    "notes": definition["notes"],
                }
            )
            continue

        counts = {"case": 0, "control": 0, "missing": 0}
        for eid in usable_ids:
            row = phenotypes.get(eid)
            value = "" if row is None else row.get(phenotype_column, "")
            status = case_status(value, positive_values)
            counts[status] += 1
            group_subject_status[disease_group][eid].append(status)

            if status == "case":
                subject_cases.append(
                    {
                        "cohort_id": cohort.cohort_id,
                        "eid": eid,
                        "disease_group": disease_group,
                        "disease_id": disease_id,
                        "disease_label": definition["disease_label"],
                        "phenotype_column": phenotype_column,
                        "phenotype_value": value,
                    }
                )

        nonmissing_n = counts["case"] + counts["control"]
        prevalence = counts["case"] / nonmissing_n if nonmissing_n else math.nan
        disease_counts.append(
            {
                "cohort_id": cohort.cohort_id,
                "cohort_label": cohort.cohort_label,
                "analysis_role": cohort.analysis_role,
                "disease_group": disease_group,
                "disease_id": disease_id,
                "disease_label": definition["disease_label"],
                "phenotype_column": phenotype_column,
                "usable_total_n": len(usable_ids),
                "case_n": counts["case"],
                "control_n": counts["control"],
                "missing_n": counts["missing"],
                "nonmissing_n": nonmissing_n,
                "prevalence": "" if math.isnan(prevalence) else f"{prevalence:.6g}",
                "ready": "true",
                "issue": "",
                "notes": definition["notes"],
            }
        )

    group_counts: List[Dict[str, object]] = []
    for disease_group, eid_statuses in sorted(group_subject_status.items()):
        any_case_n = 0
        no_case_n = 0
        missing_all_n = 0
        for eid in usable_ids:
            statuses = eid_statuses.get(eid, [])
            if statuses and any(status == "case" for status in statuses):
                any_case_n += 1
            elif statuses and any(status != "missing" for status in statuses):
                no_case_n += 1
            else:
                missing_all_n += 1
        used_columns = [
            row["phenotype_column"]
            for row in disease_counts
            if row["disease_group"] == disease_group and row["phenotype_column"]
        ]
        group_counts.append(
            {
                "cohort_id": cohort.cohort_id,
                "cohort_label": cohort.cohort_label,
                "analysis_role": cohort.analysis_role,
                "disease_group": disease_group,
                "usable_total_n": len(usable_ids),
                "any_case_n": any_case_n,
                "no_case_n": no_case_n,
                "missing_all_n": missing_all_n,
                "disease_columns_used": ";".join(used_columns),
            }
        )

    subject_stem = {
        "gwas_core_57485": "subjects_with_disease_in_gwas_core",
        "embedding_59057": "subjects_with_disease_in_embedding_cohort",
    }.get(cohort.cohort_id, f"subjects_with_disease_{sanitize_name(cohort.cohort_id)}")

    write_tsv(
        out_dir / f"{subject_stem}.tsv",
        subject_cases,
        (
            "cohort_id",
            "eid",
            "disease_group",
            "disease_id",
            "disease_label",
            "phenotype_column",
            "phenotype_value",
        ),
    )

    return {
        "cohort": cohort,
        "usable_ids": usable_ids,
        "usable_set": usable_set,
        "disease_counts": disease_counts,
        "group_counts": group_counts,
        "subject_cases": subject_cases,
    }


def sanitize_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_")


def write_summary_report(
    path: Path,
    phenotype_table: Path,
    definitions_path: Path,
    results: Sequence[Dict[str, object]],
    disease_counts: Sequence[Dict[str, object]],
    group_counts: Sequence[Dict[str, object]],
) -> None:
    lines = [
        "# rsfMRI Disease Count Summary",
        "",
        "## Inputs",
        "",
        f"- Phenotype table: `{phenotype_table}`",
        f"- Disease definitions: `{definitions_path}`",
        "",
        "## Cohorts",
        "",
    ]

    for result in results:
        cohort = result["cohort"]
        assert isinstance(cohort, Cohort)
        lines.append(
            f"- `{cohort.cohort_id}`: {cohort.cohort_label}; "
            f"role={cohort.analysis_role}; usable_n={len(result['usable_ids'])}; "
            f"list=`{cohort.usable_list}`"
        )

    lines.extend(
        [
            "",
            "## Disease Counts",
            "",
            "| cohort | group | disease | column | case_n | control_n | missing_n | prevalence | issue |",
            "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in disease_counts:
        lines.append(
            "| {cohort_id} | {disease_group} | {disease_id} | {phenotype_column} | "
            "{case_n} | {control_n} | {missing_n} | {prevalence} | {issue} |".format(**row)
        )

    lines.extend(
        [
            "",
            "## Group Counts",
            "",
            "| cohort | group | any_case_n | no_case_n | missing_all_n | disease_columns_used |",
            "| --- | --- | ---: | ---: | ---: | --- |",
        ]
    )
    for row in group_counts:
        lines.append(
            "| {cohort_id} | {disease_group} | {any_case_n} | {no_case_n} | "
            "{missing_all_n} | {disease_columns_used} |".format(**row)
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = Path(args.project_root).resolve()
    phenotype_table = resolve_path(args.phenotype_table, root)
    definitions_path = resolve_path(args.definitions, root)
    out_dir = resolve_path(args.output_dir, root)
    out_dir.mkdir(parents=True, exist_ok=True)

    phenotypes, phenotype_columns = load_phenotypes(phenotype_table)
    definitions = load_definitions(definitions_path)
    cohorts = load_cohorts(args, root)

    results = [
        summarize_cohort(cohort, phenotypes, phenotype_columns, definitions, out_dir)
        for cohort in cohorts
    ]

    disease_counts = [row for result in results for row in result["disease_counts"]]
    group_counts = [row for result in results for row in result["group_counts"]]
    subject_cases = [row for result in results for row in result["subject_cases"]]
    missing_definitions = [row for row in disease_counts if row["ready"] != "true"]

    write_tsv(
        out_dir / "rsfmri_disease_counts.tsv",
        disease_counts,
        (
            "cohort_id",
            "cohort_label",
            "analysis_role",
            "disease_group",
            "disease_id",
            "disease_label",
            "phenotype_column",
            "usable_total_n",
            "case_n",
            "control_n",
            "missing_n",
            "nonmissing_n",
            "prevalence",
            "ready",
            "issue",
            "notes",
        ),
    )
    write_tsv(
        out_dir / "rsfmri_group_counts.tsv",
        group_counts,
        (
            "cohort_id",
            "cohort_label",
            "analysis_role",
            "disease_group",
            "usable_total_n",
            "any_case_n",
            "no_case_n",
            "missing_all_n",
            "disease_columns_used",
        ),
    )
    write_tsv(
        out_dir / "rsfmri_disease_missing_definitions.tsv",
        missing_definitions,
        (
            "cohort_id",
            "cohort_label",
            "analysis_role",
            "disease_group",
            "disease_id",
            "disease_label",
            "phenotype_column",
            "usable_total_n",
            "case_n",
            "control_n",
            "missing_n",
            "nonmissing_n",
            "prevalence",
            "ready",
            "issue",
            "notes",
        ),
    )
    write_tsv(
        out_dir / "subjects_with_disease_all_cohorts.tsv",
        subject_cases,
        (
            "cohort_id",
            "eid",
            "disease_group",
            "disease_id",
            "disease_label",
            "phenotype_column",
            "phenotype_value",
        ),
    )

    result_by_id = {result["cohort"].cohort_id: result for result in results}
    embedding_not_core: List[Dict[str, object]] = []
    if {"gwas_core_57485", "embedding_59057"}.issubset(result_by_id):
        core_ids = result_by_id["gwas_core_57485"]["usable_set"]
        embedding_not_core = [
            row
            for row in result_by_id["embedding_59057"]["subject_cases"]
            if row["eid"] not in core_ids
        ]
    write_tsv(
        out_dir / "disease_positive_embedding_not_gwas_core.tsv",
        embedding_not_core,
        (
            "cohort_id",
            "eid",
            "disease_group",
            "disease_id",
            "disease_label",
            "phenotype_column",
            "phenotype_value",
        ),
    )

    write_summary_report(
        out_dir / "rsfmri_disease_summary_report.md",
        phenotype_table,
        definitions_path,
        results,
        disease_counts,
        group_counts,
    )

    print("rsfMRI disease summary complete.")
    print(f"Cohorts: {len(cohorts)}")
    print(f"Phenotype participants: {len(phenotypes)}")
    print(f"Disease definitions: {len(definitions)}")
    print(f"Ready disease rows: {sum(row['ready'] == 'true' for row in disease_counts)}")
    print(f"Missing phenotype-definition rows: {len(missing_definitions)}")
    print(f"Disease-positive embedding-not-core rows: {len(embedding_not_core)}")
    print(f"Output directory: {out_dir}")


if __name__ == "__main__":
    main()
