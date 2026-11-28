# Comprehensive simulation examples

Run from the repository root:

```powershell
python .\comprehensive\generate_difficulty_set.py --out ..\simulator_runs\comprehensive_01
```

The three configurations generate one tree each. They vary family separation and late gene duplication together:

| Level | Config | Family split age range | Late duplication rate | Expected effect |
| --- | --- | ---: | ---: | --- |
| Easy | `config_random_easy.json` | 1.50–1.80 | 0.12 | More separated families; fewer extra copies. |
| Medium | `config_random_medium.json` | 1.18–1.38 | 0.35 | Intermediate separation and copy number. |
| Hard | `config_random_topology.json` | 1.04–1.12 | 0.70 | Shorter branches between families; more extra copies. |

Shared configuration fields:

| Parameter | Value | Effect |
| --- | ---: | --- |
| `k` | 4 | Number of ancestral gene families. |
| `species_count` | 12 | Number of species in the sampled species tree. |
| `species_root_age` | 1.0 | Time from the species-tree root to the present. |
| `loss_rate` | 0.02 | Gene loss event rate per lineage per unit of time. |
| `substitution_rate` | 0.3 | Baseline sequence substitutions per site per unit of time. |
| `branch_rate_multipliers` | S03 and S09: 3.0 | Multiplies sequence evolution on those species branches; changes branch lengths. |
| `sequence_length` | 250 | Nucleotides per simulated sequence. |
| `seed` | 20261126 | Starting value for reproducible random sampling. |
| `replicates` | 1 | Number of outputs per configuration in the bundle. |
| `require_all_species_families` | true | Requires every family to appear in every species. |
| `max_attempts` | 200 | Maximum candidate simulations for each output. |
| `max_events_per_attempt` | 100000 | Cap on event nodes in one candidate simulation. |
| `stem_time` | 0.2 | Time before the oldest family split. |

Each generated bundle has a `simulation_summary.md`; its subdirectories contain sequences, the true gene tree, labels, and event records.
