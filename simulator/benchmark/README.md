# Benchmark simulation

Run from the repository root with a new output path:

```powershell
python .\benchmark\generate_benchmark.py --out ..\simulator_runs\benchmark_01
```

The settings are in `benchmark_config.json`. A run produces two experiment series:

| Series | Conditions | Varying factor | Fixed factors |
| --- | --- | --- | --- |
| `difficulty` | k = 2, 3, 4 × easy, medium, hard | Family split ages and evolutionary rates on selected species branches | Late duplication rate 0.2; loss rate 0 |
| `size` | k = 4 × low, middle, high | Late duplication rate | Family split ages 1.25–1.40; species rate multipliers 1; loss rate 0 |

## Configuration parameters

| Parameter | Default | Effect |
| --- | --- | --- |
| `k_values` | 2, 3, 4 | Number of ancestral families in the difficulty series. |
| `replicates_per_cell` | 5 | Accepted trees per k × condition. |
| `seed_base` | 20261201 | Starting value for reproducible random sampling. |
| `species_count` | 6 | Number of tips in each sampled species tree. |
| `species_root_age` | 1.0 | Time from the species-tree root to the present. |
| `late_duplication_rate` | 0.2 | Gene duplication rate per lineage per unit of time in the difficulty series; affects copy number and within-family topology. |
| `min_duplicated_families_by_k` | 2: 2; 3: 2; 4: 3 | Minimum number of families with a late duplication in each difficulty sample. |
| `loss_rate` | 0 | Gene loss rate per lineage per unit of time. |
| `substitution_rate` | 0.3 | Baseline substitutions per site per unit of time. |
| `sequence_length` | 250 | Nucleotides in each simulated sequence. |
| `require_all_species_families` | true | Requires every family to appear in every species. |
| `max_attempts` | 300 | Maximum sampled candidates per accepted tree. |
| `max_events_per_attempt` | 100000 | Event-node cap for one candidate. |
| `stem_time` | 0.2 | Time before the oldest family split. |

Difficulty condition settings:

| Condition | `family_age_min/max` | `species_rate_multipliers` | Accepted sequence overlap | Accepted sequence margin |
| --- | ---: | --- | ---: | ---: |
| Easy | 1.50–1.80 | none | 0–0.03 | ≥ 0.03 |
| Medium | 1.18–1.38 | S03 × 2 | 0.02–0.65 | −0.06–0.05 |
| Hard | 1.04–1.12 | S03 × 3; S05 × 3 | 0.50–1.00 | ≤ −0.05 |

`family_age_min/max` sets the range for splits between ancestral families. Older splits generally produce greater sequence separation. `species_rate_multipliers` changes sequence branch lengths on the named species branches. The size series uses `size_experiment.k_values`, `family_age_min/max`, `species_rate_multipliers`, and `late_duplication_rates` (low 0.12, middle 0.35, high 0.70).

## Replicates and random seeds

For replicate number `r` (starting at 1), the generator assigns the starting seed as follows:

| Series | Starting seed |
| --- | --- |
| Difficulty | `seed_base + 10000 × k + (r − 1)` |
| Size | `seed_base + 10000 × k + (r − 1) + 1000000` |

Different replicate numbers use different starting seeds. Within a series, the same `k` and replicate number share a starting seed across conditions. This pairs the first candidate's sampled species and ancestral-family topologies. Condition settings then affect branch lengths and late duplication events. The generator continues sampling from the same random stream until a candidate meets that condition's selection criteria, so the final accepted gene trees can have different topologies or accepted attempt numbers. The `seed` and `accepted_attempt` columns in `manifest.tsv` record these details. Reusing the same settings and seed reproduces the same batch.

## Measurements and output

The generator compares every sequence pair using raw p-distance, the fraction of differing nucleotide positions. It also records the corresponding distances on the true tree.

| Measure | Definition | Role |
| --- | --- | --- |
| `sequence_overlap_fraction` | Fraction of between-family pairs at or below the 90th percentile of within-family sequence distances. | Difficulty selection and reported score. |
| `sequence_margin` | 10th percentile of between-family distances minus 90th percentile of within-family distances. | Difficulty selection; positive values indicate greater separation. |
| `duplicated_families` | Number of families containing at least one late duplication. | Topology-diversity selection in the difficulty series. |

The **Accepted sequence margin** column in the difficulty table gives the allowed range for `sequence_margin`. For example, if the 90th percentile of within-family p-distances is 0.40 and the 10th percentile of between-family p-distances is 0.46, the margin is `0.46 − 0.40 = +0.06`. A positive margin means these two percentile cutoffs are separated; a negative margin means the between-family cutoff is below the within-family 90th percentile. These are distance differences, measured as fractions of sequence positions.

| Output | Contents |
| --- | --- |
| `benchmark_summary.md` | Mean measures for each condition. |
| `manifest.tsv` | Per-tree measures, accepted attempt, seed, and output directory. |
| `rep_*/sequences.fasta` | Sequences for tree inference and clustering. |
| `rep_*/true_gene_tree.nwk` | Ground-truth gene tree with branch lengths. |
| `rep_*/labels.tsv` | Ground-truth family and species for each sequence. |
| `rep_*/sampled_species_tree.nwk` | Sampled six-species tree. |
| `rep_*/sampled_family_topology.nwk` | Sampled topology connecting the k ancestral families. |
| `rep_*/metadata.json` | Full configuration and tree/sequence measurements. |

For a smaller run, add `--replicates 1`. Change `seed_base` for an independent batch.
