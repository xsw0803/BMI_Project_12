#!/usr/bin/env python3
"""Generate easy, medium, and hard simulations in one portable output folder."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from random_topology_simulator import check_config, read_ancestor, save_one


LEVELS = (
    ("easy", "config_random_easy.json"),
    ("medium", "config_random_medium.json"),
    ("hard", "config_random_topology.json"),
)


def summary_md(records: list[tuple[str, str, dict, dict]]) -> str:
    first = records[0][2]
    lines = [
        "# Simulation summary",
        "",
        "| Setting | Value |",
        "| --- | ---: |",
        "| Replicates per level | 1 |",
        f"| Species | {first['species_count']} |",
        f"| Families (k) | {first['k']} |",
        f"| Seed | {first['seed']} |",
        f"| Sequence length | {records[0][3]['sequence_length']} nt |",
        f"| Species root age | {first['species_root_age']} |",
        f"| Substitution rate | {first['substitution_rate']} |",
        f"| Loss rate | {first['loss_rate']} |",
        "",
        "| Level | Directory | Family split age range | Late duplication rate | Leaves | k | Late duplications | Losses | Accepted attempt |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for level, folder, cfg, meta in records:
        lines.append(
            f"| {level} | [{folder}]({folder}/) | "
            f"{cfg['family_age_min']}–{cfg['family_age_max']} | "
            f"{cfg['late_duplication_rate']} | "
            f"{meta['leaf_count']} | {meta['family_count']} | "
            f"{meta['late_duplication_events']} | {meta['loss_events']} | "
            f"{meta['accepted_attempt']} |"
        )
    lines.extend(("", "Each directory contains simulated sequences in `sequences.fasta`, with evaluation truth in `true_gene_tree.nwk` and `labels.tsv`.", ""))
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True,
                        help="New folder that will contain all three difficulty outputs.")
    parser.add_argument("--config-dir", type=Path, default=Path(__file__).resolve().parent,
                        help="Directory containing the three config JSON files.")
    args = parser.parse_args()
    output = args.out.resolve()
    if output.exists():
        parser.error(f"Output folder already exists: {output}")

    inputs = []
    for level, filename in LEVELS:
        path = args.config_dir / filename
        cfg = json.loads(path.read_text(encoding="utf-8"))
        check_config(cfg)
        if cfg["replicates"] != 1:
            parser.error(f"{filename}: replicates must be 1 for this bundle")
        ancestor = read_ancestor(cfg, path)
        inputs.append((level, path, cfg, ancestor))

    shared_keys = (
        "k", "species_count", "species_root_age", "substitution_rate",
        "loss_rate", "branch_rate_multipliers", "sequence_length",
        "require_all_species_families",
    )
    baseline = inputs[0][2]
    for level, _, cfg, ancestor in inputs[1:]:
        for key in shared_keys:
            if cfg[key] != baseline[key]:
                parser.error(f"{level}: {key} differs from easy; update the summary design")
        if len(ancestor or "") != len(inputs[0][3] or ""):
            parser.error(f"{level}: ancestral sequence length differs from easy")

    records = []
    for level, path, cfg, ancestor in inputs:
        folder_name = f"simulated_{level}_01"
        folder = output / folder_name
        meta = save_one(cfg, folder, 0, ancestor)
        portable_cfg = dict(cfg)
        if "ancestor_fasta" in portable_cfg:
            portable_cfg.pop("ancestor_fasta")
            portable_cfg["ancestor_sequence"] = ancestor
        (folder / "config_used.json").write_text(
            json.dumps(portable_cfg, indent=2) + "\n", encoding="utf-8"
        )
        records.append((level, folder_name, cfg, meta))
        print(f"{folder_name}: {meta['leaf_count']} leaves, "
              f"{meta['family_count']} families, "
              f"{meta['late_duplication_events']} late duplications, "
              f"{meta['loss_events']} losses; attempt {meta['accepted_attempt']}")

    (output / "simulation_summary.md").write_text(summary_md(records), encoding="utf-8")
    print(f"Summary: {output / 'simulation_summary.md'}")


if __name__ == "__main__":
    main()
