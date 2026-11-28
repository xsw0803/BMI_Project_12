"""Shared gene-lineage simulation and JC69 sequence evolution."""

from __future__ import annotations

import math
import random
from collections import Counter
from dataclasses import dataclass, field
from typing import Any


@dataclass

class Species:
    name: str
    incoming_time: float = 0.0
    children: list[Species] = field(default_factory=list)
    age: float = 0.0


@dataclass

class Gene:
    node_id: int
    age: float
    event: str
    family: str | None = None
    species: str | None = None
    length: float = 0.0  # incoming edge, substitutions/site
    children: list[Gene] = field(default_factory=list)
    label: str | None = None


class AttemptOverflow(Exception):
    """A high birth rate generated too many event nodes in one attempt."""


class Simulator:
    def __init__(self, config: dict[str, Any], species: Species, rng: random.Random):
        self.cfg = config
        self.species = species
        self.rng = rng
        self.next_id = 0
        self.events: list[dict[str, Any]] = []
        self.max_events = config.get("max_events_per_attempt", 100000)

    def node(self, age: float, event: str, parent: Gene | None,
             family: str | None = None, species: str | None = None,
             multiplier: float = 1.0) -> Gene:
        self.next_id += 1
        if self.next_id > self.max_events:
            raise AttemptOverflow("Event count exceeded max_events_per_attempt")
        duration = 0.0 if parent is None else parent.age - age
        if duration < -1e-9:
            raise RuntimeError("Event times are in the wrong order.")
        length = max(0.0, duration) * self.cfg["substitution_rate"] * multiplier
        item = Gene(self.next_id, age, event, family, species, length)
        if parent is not None:
            parent.children.append(item)
        if event in ("early_duplication", "late_duplication", "speciation"):
            self.events.append({"event": event, "node_id": item.node_id,
                                "age": age, "family": family, "species_branch": species})
        return item

    def factor(self, family: str, species_branch: str | None = None) -> float:
        settings = self.cfg["branch_rate_multipliers"]
        family_factor = settings.get("family", {}).get(family, 1.0)
        species_factor = settings.get("species", {}).get(species_branch, 1.0)
        return family_factor * species_factor

    def family_at_species_root(self, parent: Gene, family: str) -> None:
        self.evolve_species(parent, self.species, family, None)

    def evolve_species(self, parent: Gene, species_node: Species, family: str,
                       branch_name: str | None) -> None:
        """One lineage from parent.age down an incoming species branch."""
        birth = self.cfg["late_duplication_rate"] if branch_name is not None else 0.0
        death = self.cfg["loss_rate"] if branch_name is not None else 0.0
        rate = birth + death
        if rate:
            wait = self.rng.expovariate(rate)
            event_age = parent.age - wait
            if event_age > species_node.age + 1e-12:
                if self.rng.random() < death / rate:
                    self.events.append({"event": "loss", "age": event_age,
                                        "family": family, "species_branch": branch_name,
                                        "parent_node_id": parent.node_id})
                    return
                duplicate = self.node(event_age, "late_duplication", parent,
                                      family, branch_name, self.factor(family, branch_name))
                self.evolve_species(duplicate, species_node, family, branch_name)
                self.evolve_species(duplicate, species_node, family, branch_name)
                return
        if not species_node.children:
            self.node(species_node.age, "tip", parent, family, species_node.name,
                      self.factor(family, branch_name))
            return
        split = self.node(species_node.age, "speciation", parent, family,
                          species_node.name, self.factor(family, branch_name))
        for child in species_node.children:
            self.evolve_species(split, child, family, child.name)


def prune(node: Gene) -> Gene | None:
    if node.event == "tip":
        return node
    node.children = [survivor for child in node.children
                     if (survivor := prune(child)) is not None]
    if not node.children:
        return None
    if len(node.children) == 1:
        only = node.children[0]
        only.length += node.length
        return only
    return node


def tips(root: Gene) -> list[Gene]:
    if not root.children:
        return [root]
    return [tip for child in root.children for tip in tips(child)]


def label_tips(root: Gene) -> list[Gene]:
    leaves = tips(root)
    counts: Counter[tuple[str, str]] = Counter()
    for leaf in leaves:
        assert leaf.family is not None and leaf.species is not None
        key = leaf.family, leaf.species
        counts[key] += 1
        leaf.label = f"{leaf.family}__{leaf.species}__copy{counts[key]}"
    return leaves


def newick(root: Gene) -> str:
    def render(node: Gene, is_root: bool) -> str:
        if node.children:
            label = "(" + ",".join(render(child, False) for child in node.children) + ")"
        else:
            assert node.label is not None
            label = node.label
        return label if is_root else f"{label}:{node.length:.12g}"
    return render(root, True) + ";\n"


def evolve_sequence(seq: str, length: float, rng: random.Random) -> str:
    probability_change = -0.75 * math.expm1(-4.0 * length / 3.0)
    alphabet = "ACGT"
    result: list[str] = []
    for base in seq:
        if rng.random() < probability_change:
            result.append(alphabet[(alphabet.index(base) + rng.randrange(1, 4)) % 4])
        else:
            result.append(base)
    return "".join(result)
