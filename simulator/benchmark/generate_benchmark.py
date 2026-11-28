#!/usr/bin/env python3
"""Generate a k × difficulty benchmark with automatic distance-overlap checks."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from statistics import mean
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from random_topology_simulator import check_config, save_one  # noqa: E402


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower) if lower != upper else ordered[lower]


def pair_metrics(within: list[float], between: list[float], prefix: str) -> dict[str, Any]:
    if not within or not between:
        raise ValueError("Distance metric needs both within- and between-family pairs")
    within_q90 = quantile(within, 0.9)
    between_q10 = quantile(between, 0.1)
    overlap_fraction = sum(value <= within_q90 for value in between) / len(between)
    return {
        f"{prefix}_within_q90": round(within_q90, 8),
        f"{prefix}_between_q10": round(between_q10, 8),
        f"{prefix}_margin": round(between_q10 - within_q90, 8),
        f"{prefix}_overlap_fraction": round(overlap_fraction, 6),
    }


def distance_diagnostics(root: Any, leaves: list[Any],
                         sequences: dict[str, str]) -> dict[str, Any]:
    within: list[float] = []
    between: list[float] = []

    def visit(node: Any, depth: float) -> list[tuple[str, float]]:
        if not node.children:
            return [(node.family, depth)]
        child_groups = [visit(child, depth + child.length) for child in node.children]
        for index, group in enumerate(child_groups):
            for other in child_groups[index + 1:]:
                for family_a, depth_a in group:
                    for family_b, depth_b in other:
                        distance = depth_a + depth_b - 2 * depth
                        (within if family_a == family_b else between).append(distance)
        return [leaf for group in child_groups for leaf in group]

    visit(root, 0.0)
    sequence_within: list[float] = []
    sequence_between: list[float] = []
    for index, first in enumerate(leaves):
        for second in leaves[index + 1:]:
            left, right = sequences[first.label], sequences[second.label]
            distance = sum(a != b for a, b in zip(left, right)) / len(left)
            (sequence_within if first.family == second.family else sequence_between).append(distance)
    return {
        "within_pairs": len(within),
        "between_pairs": len(between),
        **pair_metrics(within, between, "tree"),
        **pair_metrics(sequence_within, sequence_between, "sequence"),
    }


def duplication_counts(root: Any) -> dict[str, int]:
    counts: dict[str, int] = {}

    def visit(node: Any) -> None:
        if node.event == "late_duplication":
            counts[node.family] = counts.get(node.family, 0) + 1
        for child in node.children:
            visit(child)

    visit(root)
    return counts


def make_sim_config(settings: dict[str, Any], k: int,
                    condition: dict[str, Any], seed: int,
                    late_duplication_rate: float | None = None) -> dict[str, Any]:
    cfg = {key: settings[key] for key in (
        "species_count", "species_root_age", "late_duplication_rate",
        "loss_rate", "substitution_rate", "sequence_length",
        "require_all_species_families", "max_attempts",
        "max_events_per_attempt", "stem_time",
    )}
    cfg.update({
        "k": k,
        "family_age_min": condition["family_age_min"],
        "family_age_max": condition["family_age_max"],
        "branch_rate_multipliers": {
            "species": condition["species_rate_multipliers"],
            "family": {},
        },
        "seed": seed,
        "replicates": 1,
    })
    if late_duplication_rate is not None:
        cfg["late_duplication_rate"] = late_duplication_rate
    check_config(cfg)
    return cfg


def validate_settings(settings: dict[str, Any]) -> None:
    k_values = settings["k_values"]
    if not isinstance(k_values, list) or not k_values or any(type(k) is not int or k < 2 for k in k_values):
        raise ValueError("k_values must be a nonempty list of integers >= 2")
    if len(set(k_values)) != len(k_values):
        raise ValueError("k_values must be unique")
    minimums = settings["min_duplicated_families_by_k"]
    if set(minimums) != {str(k) for k in k_values} or any(
        type(minimums[str(k)]) is not int or not 0 <= minimums[str(k)] <= k
        for k in k_values
    ):
        raise ValueError("min_duplicated_families_by_k must give an integer from 0 to k for each k")
    if type(settings["replicates_per_cell"]) is not int or settings["replicates_per_cell"] <= 0:
        raise ValueError("replicates_per_cell must be positive")
    if type(settings["seed_base"]) is not int:
        raise ValueError("seed_base must be an integer")
    conditions = settings["conditions"]
    if not isinstance(conditions, dict) or set(conditions) != {"easy", "medium", "hard"}:
        raise ValueError("conditions must contain easy, medium, hard")
    for level, condition in conditions.items():
        low, high = condition["overlap_fraction_min"], condition["overlap_fraction_max"]
        if not 0 <= low <= high <= 1:
            raise ValueError(f"{level}: overlap fraction range must lie within [0,1]")
        if condition.get("distance_margin_min", -math.inf) > condition.get("distance_margin_max", math.inf):
            raise ValueError(f"{level}: distance margin minimum exceeds maximum")
        make_sim_config(settings, k_values[0], condition, settings["seed_base"])
    size = settings["size_experiment"]
    if not isinstance(size["k_values"], list) or not size["k_values"]:
        raise ValueError("size_experiment.k_values must be nonempty")
    if any(type(k) is not int or k < 2 for k in size["k_values"]):
        raise ValueError("size_experiment.k_values must contain integers >= 2")
    if len(set(size["k_values"])) != len(size["k_values"]):
        raise ValueError("size_experiment.k_values must be unique")
    rates = size["late_duplication_rates"]
    if not isinstance(rates, dict) or set(rates) != {"low", "middle", "high"}:
        raise ValueError("size_experiment.late_duplication_rates needs low, middle, high")
    for level, rate in rates.items():
        make_sim_config(settings, size["k_values"][0], size, settings["seed_base"], rate)


def summary(rows: list[dict[str, Any]], settings: dict[str, Any]) -> str:
    lines = [
        "# Benchmark summary", "",
        "| Setting | Value |", "| --- | ---: |",
        f"| Replicates per condition | {settings['replicates_per_cell']} |",
        f"| Species | {settings['species_count']} |",
        f"| Sequence length | {settings['sequence_length']} nt |",
        f"| Difficulty-series late duplication rate | {settings['late_duplication_rate']} |",
        f"| Loss rate | {settings['loss_rate']} |", "",
        "Sequence overlap fraction = proportion of between-family leaf pairs whose raw sequence p-distance",
        "is no greater than the 90th percentile of within-family distances.", "",
        "## Difficulty series", "",
        "Family split ages and selected species rates vary. The late duplication rate is fixed.",
        "Selection criteria: sequence overlap, sequence margin, and duplicated-family count.",
        "", "| k | Difficulty | Replicates | Mean leaves | Mean overlap | Mean margin | Mean late duplications | Mean duplicated families |",
        "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for k in settings["k_values"]:
        for level in ("easy", "medium", "hard"):
            cell = [row for row in rows if row["experiment"] == "difficulty" and row["k"] == k and row["condition"] == level]
            lines.append(
                f"| {k} | {level} | {len(cell)} | "
                f"{mean(row['leaves'] for row in cell):.1f} | "
                f"{mean(row['sequence_overlap_fraction'] for row in cell):.3f} | "
                f"{mean(row['sequence_margin'] for row in cell):.3f} | "
                f"{mean(row['late_duplications'] for row in cell):.1f} | "
                f"{mean(row['duplicated_families'] for row in cell):.1f} |"
            )
    lines.extend(("", "## Size series", "",
                  "Family split ages and species rates are fixed. The late duplication rate varies.",
                  "", "| k | Duplication setting | Rate | Replicates | Mean leaves | Mean late duplications |",
                  "| ---: | --- | ---: | ---: | ---: | ---: |"))
    for k in settings["size_experiment"]["k_values"]:
        for level, rate in settings["size_experiment"]["late_duplication_rates"].items():
            cell = [row for row in rows if row["experiment"] == "size" and row["k"] == k and row["condition"] == level]
            lines.append(f"| {k} | {level} | {rate} | {len(cell)} | "
                         f"{mean(row['leaves'] for row in cell):.1f} | "
                         f"{mean(row['late_duplications'] for row in cell):.1f} |")
    lines.extend(("", "Per-tree values and accepted seeds are in `manifest.tsv`.",
                  "`true_gene_tree.nwk` and `labels.tsv` are evaluation truth.", ""))
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("benchmark_config.json"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--replicates", type=int, help="Override replicates_per_cell for a pilot run")
    args = parser.parse_args()
    settings = json.loads(args.config.read_text(encoding="utf-8"))
    if args.replicates is not None:
        settings["replicates_per_cell"] = args.replicates
    validate_settings(settings)
    output = args.out.resolve()
    if output.exists():
        parser.error(f"Output directory already exists: {output}")
    output.mkdir(parents=True)
    (output / "benchmark_config_used.json").write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")

    rows: list[dict[str, Any]] = []
    def generate_cell(experiment: str, k: int, level: str,
                      condition: dict[str, Any], rate: float,
                      seed_offset: int) -> None:
        for rep in range(settings["replicates_per_cell"]):
            seed = settings["seed_base"] + 10000 * k + rep + seed_offset
            cfg = make_sim_config(settings, k, condition, seed, rate)
            folder = output / experiment / f"k{k}" / level / f"rep_{rep + 1:03d}"

            def accept(root: Any, leaves: list[Any], sequences: dict[str, str]) -> tuple[bool, dict[str, Any]]:
                diagnostics = distance_diagnostics(root, leaves, sequences)
                counts = duplication_counts(root)
                diagnostics["late_duplications_by_family"] = counts
                diagnostics["duplicated_families"] = len(counts)
                score = diagnostics["sequence_overlap_fraction"]
                margin = diagnostics["sequence_margin"]
                keep = (condition.get("overlap_fraction_min", 0) <= score <= condition.get("overlap_fraction_max", 1)
                        and condition.get("distance_margin_min", -math.inf) <= margin
                        <= condition.get("distance_margin_max", math.inf))
                if experiment == "difficulty":
                    keep = keep and len(counts) >= settings["min_duplicated_families_by_k"][str(k)]
                return keep, diagnostics

            try:
                meta = save_one(cfg, folder, 0, None, sample_filter=accept)
            except ValueError as exc:
                raise ValueError(f"{experiment}, k={k}, {level}, rep={rep + 1}: {exc}") from exc
            row = {
                "experiment": experiment, "k": k, "condition": level,
                "replicate": rep + 1, "seed": meta["seed"],
                "accepted_attempt": meta["accepted_attempt"],
                "late_duplication_rate": rate,
                "leaves": meta["leaf_count"],
                "late_duplications": meta["late_duplication_events"],
                "duplicated_families": meta["distance_diagnostics"]["duplicated_families"],
                "losses": meta["loss_events"],
                **{key: value for key, value in meta["distance_diagnostics"].items()
                   if key != "late_duplications_by_family"},
                "late_duplications_by_family": json.dumps(
                    meta["distance_diagnostics"]["late_duplications_by_family"], sort_keys=True
                ),
                "folder": folder.relative_to(output).as_posix(),
            }
            rows.append(row)
            print(f"{experiment} k={k} {level} rep_{rep + 1:03d}: "
                  f"{row['leaves']} leaves, overlap={row['sequence_overlap_fraction']:.3f}, "
                  f"margin={row['sequence_margin']:.3f}, "
                  f"attempt={row['accepted_attempt']}")

    for k in settings["k_values"]:
        for level in ("easy", "medium", "hard"):
            generate_cell("difficulty", k, level, settings["conditions"][level],
                          settings["late_duplication_rate"], 0)
    size = settings["size_experiment"]
    for k in size["k_values"]:
        for level, rate in size["late_duplication_rates"].items():
            generate_cell("size", k, level, size, rate, 1000000)

    with (output / "manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (output / "benchmark_summary.md").write_text(summary(rows, settings), encoding="utf-8")
    print(f"Finished {len(rows)} trees in {output}")


if __name__ == "__main__":
    main()
