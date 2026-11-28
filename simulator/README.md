# Gene family simulator

This is a portable Python 3.10+ simulator using only the standard library. It samples species and ancestral-family topologies, simulates gene duplication and loss, and evolves nucleotide sequences on the resulting gene tree.

| Path | Purpose |
| --- | --- |
| `engine_core.py` | Gene-lineage events and JC69 sequence evolution. |
| `random_topology_simulator.py` | Random topology sampling and per-tree output. |
| [`comprehensive/`](comprehensive/) | Three example levels that vary family separation and duplication together. |
| [`benchmark/`](benchmark/) | Controlled difficulty and tree-size experiments with measured sequence distances. |

Run either generator from this directory. The examples write output to a sibling `simulator_runs` directory:

```powershell
python .\comprehensive\generate_difficulty_set.py --out ..\simulator_runs\comprehensive_01
python .\benchmark\generate_benchmark.py --out ..\simulator_runs\benchmark_01
```

The generated `sequences.fasta` files are the sequence data. `true_gene_tree.nwk` and `labels.tsv` record the evaluation truth. See each workflow's README for parameters and output details.

## Sequence evolution model

Nucleotide substitutions along each gene-tree branch follow the Jukes–Cantor 1969 (JC69) model, which assigns equal frequencies to A, C, G, and T and equal rates to all substitutions between them. The implementation is in `engine_core.py`.

**Reference:** Jukes, T. H., & Cantor, C. R. (1969). Evolution of protein molecules. In H. N. Munro (Ed.), *Mammalian Protein Metabolism* (Vol. 3, pp. 21–132). Academic Press. [Original chapter PDF](https://evoluscope.fr/phylographe/biblio/JukesCantor1969.pdf).
