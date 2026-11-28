#!/usr/bin/env python3
"""Sample species and ancestral-family topologies and evolve gene copies.

Writes nucleotide sequences, the true gene tree, and family labels for evaluation.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from engine_core import (
    AttemptOverflow, Gene, Simulator, Species, evolve_sequence, label_tips,
    newick, prune,
)


EXAMPLE = {
    "k": 4,
    "species_count": 12,
    "species_root_age": 1.0,
    "family_age_min": 1.04,
    "family_age_max": 1.3,
    "late_duplication_rate": 0.35,
    "loss_rate": 0.02,
    "substitution_rate": 0.3,
    "branch_rate_multipliers": {"species": {}, "family": {}},
    "sequence_length": 600,
    "seed": 20261125,
    "replicates": 1,
    "require_all_species_families": True,
    "max_attempts": 200,
    "max_events_per_attempt": 100000,
    "stem_time": 0.2,
}


@dataclass
class FamilyBranch:
    age: float
    family: str | None = None
    children: tuple[FamilyBranch, FamilyBranch] | None = None


def check_config(cfg: dict[str, Any]) -> None:
    for key in ("k", "species_count", "seed", "replicates", "max_attempts",
                "max_events_per_attempt"):
        if type(cfg.get(key)) is not int or (key != "seed" and cfg[key] <= 0):
            raise ValueError(f"{key} must be an integer (positive except seed)")
    if cfg["k"] < 2 or cfg["species_count"] < 2:
        raise ValueError("k and species_count must each be at least 2")
    for key in ("species_root_age", "family_age_min", "family_age_max",
                "late_duplication_rate", "loss_rate", "substitution_rate", "stem_time"):
        value = cfg.get(key)
        if type(value) not in (int, float) or not 0 <= value < float("inf"):
            raise ValueError(f"{key} must be a finite nonnegative number")
    if cfg["species_root_age"] <= 0 or cfg["substitution_rate"] <= 0:
        raise ValueError("species_root_age and substitution_rate must be positive")
    if not cfg["species_root_age"] < cfg["family_age_min"] < cfg["family_age_max"]:
        raise ValueError("family_age_min/max must be above species_root_age")
    if type(cfg.get("require_all_species_families")) is not bool:
        raise ValueError("require_all_species_families must be Boolean")
    if type(cfg.get("sequence_length")) is not int or cfg["sequence_length"] <= 0:
        raise ValueError("sequence_length must be a positive integer")
    multipliers = cfg.get("branch_rate_multipliers")
    if not isinstance(multipliers, dict) or set(multipliers) != {"species", "family"}:
        raise ValueError("branch_rate_multipliers needs species and family maps")
    allowed = {
        "species": {f"S{i:02d}" for i in range(1, cfg["species_count"] + 1)},
        "family": {f"F{i}" for i in range(1, cfg["k"] + 1)},
    }
    for kind in ("species", "family"):
        mapping = multipliers[kind]
        if not isinstance(mapping, dict):
            raise ValueError(f"{kind} multipliers must be a map")
        for name, value in mapping.items():
            if name not in allowed[kind] or type(value) not in (int, float) or not 0 < value < float("inf"):
                raise ValueError(f"invalid {kind} multiplier: {name}={value}")


def random_species_tree(count: int, root_age: float, rng: random.Random) -> Species:
    pool = [Species(f"S{i:02d}", age=0.0) for i in range(1, count + 1)]
    ages = sorted(rng.uniform(0.01 * root_age, 0.99 * root_age)
                  for _ in range(count - 2)) + [root_age]
    for number, age in enumerate(ages, 1):
        first, second = rng.sample(pool, 2)
        pool.remove(first)
        pool.remove(second)
        first.incoming_time = age - first.age
        second.incoming_time = age - second.age
        pool.append(Species(f"N{number:02d}", children=[first, second], age=age))
    return pool[0]


def species_newick(root: Species) -> str:
    def render(node: Species, is_root: bool) -> str:
        subtree = ("(" + ",".join(render(child, False) for child in node.children) + ")") if node.children else ""
        return subtree + node.name + ("" if is_root else f":{node.incoming_time:.12g}")
    return render(root, True) + ";\n"


def random_family_tree(k: int, minimum: float, maximum: float,
                       rng: random.Random) -> FamilyBranch:
    pool = [FamilyBranch(0.0, f"F{i}") for i in range(1, k + 1)]
    ages = sorted(rng.uniform(minimum, maximum) for _ in range(k - 1))
    for age in ages:
        first, second = rng.sample(pool, 2)
        pool.remove(first)
        pool.remove(second)
        pool.append(FamilyBranch(age, children=(first, second)))
    return pool[0]


def family_newick(root: FamilyBranch) -> str:
    if root.family is not None:
        return root.family
    assert root.children is not None
    return "(" + ",".join(family_newick(child) for child in root.children) + ")"


def build_gene_tree(sim: Simulator, family_tree: FamilyBranch,
                    stem_time: float) -> Gene:
    root = sim.node(family_tree.age + stem_time, "origin", None)

    def descend(branch: FamilyBranch, parent: Gene) -> None:
        if branch.family is not None:
            sim.family_at_species_root(parent, branch.family)
            return
        node = sim.node(branch.age, "early_duplication", parent)
        assert branch.children is not None
        for child in branch.children:
            descend(child, node)

    descend(family_tree, root)
    result = prune(root)
    if result is None:
        raise RuntimeError("All lineages were lost")
    result.length = 0.0
    return result


def read_ancestor(cfg: dict[str, Any], config_path: Path) -> str | None:
    if "ancestor_sequence" in cfg and "ancestor_fasta" in cfg:
        raise ValueError("Use only one of ancestor_sequence and ancestor_fasta")
    if "ancestor_fasta" in cfg:
        path = (config_path.parent / cfg["ancestor_fasta"]).resolve()
        records = [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines()]
        headers = [i for i, line in enumerate(records) if line.startswith(">")]
        if len(headers) != 1 or headers[0] != 0:
            raise ValueError("ancestor_fasta must contain exactly one FASTA record")
        sequence = "".join(records[1:]).upper()
    elif "ancestor_sequence" in cfg:
        sequence = str(cfg["ancestor_sequence"]).strip().upper()
    else:
        return None
    if not sequence or any(base not in "ACGT" for base in sequence):
        raise ValueError("Ancestor sequence must contain only A, C, G, T")
    return sequence


def sequences_from_ancestor(root: Gene, ancestor: str,
                            rng: random.Random) -> dict[str, str]:
    sequences: dict[str, str] = {}

    def visit(node: Gene, sequence: str) -> None:
        if node.label is not None:
            sequences[node.label] = sequence
        for child in node.children:
            visit(child, evolve_sequence(sequence, child.length, rng))

    visit(root, ancestor)
    return sequences


def save_one(cfg: dict[str, Any], folder: Path, replicate_index: int,
             supplied_ancestor: str | None,
             sample_filter: Callable[[Gene, list[Gene], dict[str, str]], tuple[bool, dict[str, Any]]] | None = None,
             ) -> dict[str, Any]:
    rng = random.Random(cfg["seed"] + replicate_index)
    names = [f"S{i:02d}" for i in range(1, cfg["species_count"] + 1)]
    for attempt in range(1, cfg["max_attempts"] + 1):
        species = random_species_tree(cfg["species_count"], cfg["species_root_age"], rng)
        family_tree = random_family_tree(cfg["k"], cfg["family_age_min"],
                                         cfg["family_age_max"], rng)
        sim = Simulator(cfg, species, rng)
        try:
            root = build_gene_tree(sim, family_tree, cfg["stem_time"])
        except (AttemptOverflow, RuntimeError, RecursionError):
            continue
        leaves = label_tips(root)
        observed = {(leaf.family, leaf.species) for leaf in leaves}
        if cfg["require_all_species_families"]:
            required = {(f"F{i}", name) for i in range(1, cfg["k"] + 1) for name in names}
            accepted = required <= observed
        else:
            accepted = all(any(leaf.family == f"F{i}" for leaf in leaves)
                           for i in range(1, cfg["k"] + 1))
        if accepted and len(leaves) >= 2:
            diagnostics: dict[str, Any] = {}
            if sample_filter is not None:
                ancestor = supplied_ancestor or "".join(rng.choices("ACGT", k=cfg["sequence_length"]))
                sequences = sequences_from_ancestor(root, ancestor, rng)
                keep, diagnostics = sample_filter(root, leaves, sequences)
                if not keep:
                    continue
            break
    else:
        raise ValueError("No valid sample within max_attempts")

    if sample_filter is None:
        ancestor = supplied_ancestor or "".join(rng.choices("ACGT", k=cfg["sequence_length"]))
        sequences = sequences_from_ancestor(root, ancestor, rng)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "true_gene_tree.nwk").write_text(newick(root), encoding="utf-8")
    (folder / "sampled_species_tree.nwk").write_text(species_newick(species), encoding="utf-8")
    (folder / "sampled_family_topology.nwk").write_text(
        family_newick(family_tree) + ";\n", encoding="utf-8",
    )
    (folder / "ancestor.fasta").write_text(
        ">ancestor\n" + "\n".join(ancestor[i:i + 80] for i in range(0, len(ancestor), 80)) + "\n",
        encoding="utf-8",
    )
    with (folder / "sequences.fasta").open("w", encoding="utf-8") as out:
        for leaf in leaves:
            seq = sequences[leaf.label]
            out.write(f">{leaf.label}\n")
            for start in range(0, len(seq), 80):
                out.write(seq[start:start + 80] + "\n")
    (folder / "labels.tsv").write_text(
        "leaf\tfamily\tspecies\tcopy\n" +
        "".join(f"{leaf.label}\t{leaf.family}\t{leaf.species}\t"
                f"{leaf.label.rsplit('copy', 1)[1]}\n" for leaf in leaves),
        encoding="utf-8",
    )
    (folder / "events.jsonl").write_text(
        "".join(json.dumps(event, sort_keys=True) + "\n" for event in sim.events),
        encoding="utf-8",
    )
    counts = Counter(leaf.family for leaf in leaves)
    meta = {
        "replicate": replicate_index + 1,
        "seed": cfg["seed"] + replicate_index,
        "accepted_attempt": attempt,
        "family_count": cfg["k"],
        "leaf_count": len(leaves),
        "family_sizes": dict(counts),
        "late_duplication_events": sum(e["event"] == "late_duplication" for e in sim.events),
        "loss_events": sum(e["event"] == "loss" for e in sim.events),
        "species_tree": species_newick(species).strip(),
        "family_topology": family_newick(family_tree) + ";",
        "family_split_ages": sorted(
            (e["age"] for e in sim.events if e["event"] == "early_duplication"),
            reverse=True,
        ),
        "ancestor_sequence_source": "provided" if supplied_ancestor else "random",
        "sequence_length": len(ancestor),
        "config": cfg,
    }
    if sample_filter is not None:
        meta["distance_diagnostics"] = diagnostics
    (folder / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return meta


def main() -> None:
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write-example", type=Path)
    cli.add_argument("--config", type=Path)
    cli.add_argument("--out", type=Path)
    args = cli.parse_args()
    if args.write_example:
        args.write_example.write_text(json.dumps(EXAMPLE, indent=2) + "\n", encoding="utf-8")
        print(f"Wrote {args.write_example}")
        return
    if not args.config or not args.out:
        cli.error("provide --config and --out, or --write-example")
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    check_config(cfg)
    ancestor = read_ancestor(cfg, args.config)
    for index in range(cfg["replicates"]):
        folder = args.out if cfg["replicates"] == 1 else args.out / f"rep_{index + 1:03d}"
        meta = save_one(cfg, folder, index, ancestor)
        print(f"{folder}: {meta['leaf_count']} leaves, {meta['family_count']} families, "
              f"{meta['late_duplication_events']} late duplications, "
              f"{meta['loss_events']} losses; attempt {meta['accepted_attempt']}")


if __name__ == "__main__":
    main()
