# tskit_multichrom

A Python library for efficiently storing and analyzing multiple chromosomes (contigs) using [tskit](https://tskit.dev). It allows you to work with a collection of tree sequences — one per chromosome — as a unified object while maintaining cross-chromosome sample identity.

## Why tskit_multichrom?

Population genetics simulators like **msprime** and **SLiM** typically produce separate tree sequences for each chromosome. `tskit_multichrom` provides:

- **Unified storage**: Manage multiple per-chromosome tree sequences as a single `TreesAssemblage` object
- **Cross-chromosome analysis**: Compute statistics across all chromosomes when samples are shared (e.g., diversity)
- **Efficient subsetting**: Create subsets (e.g., autosomes only) without copying underlying tree sequences
- **Format conversion**: Convert between multi-contig archives and single merged tree sequences
- **SLiM compatibility**: Seamless conversion to/from SLiM's tree sequence archives

## Quick Start

```python
import stdpopsim
import tskit_multichrom as tmc

# Simulate multiple chromosomes using stdpopsim
species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model("OutOfAfrica_2T12")
chr20 = species.get_contig("chr20", mutation_rate=model.mutation_rate)
chr21 = species.get_contig("chr21", mutation_rate=model.mutation_rate)
chr22 = species.get_contig("chr22", mutation_rate=model.mutation_rate)

samples = {'AFR': 5, 'EUR': 10}
engine = stdpopsim.get_engine("msprime")
ts20 = engine.simulate(model, chr20, samples)
ts21 = engine.simulate(model, chr21, samples)
ts22 = engine.simulate(model, chr22, samples)

# Combine into a multi-chromosome assemblage
# shared_nodes="samples" marks sample nodes as shared across chromosomes,
# which is required for cross-chromosome statistics
ta = tmc.from_tree_sequences(
    [ts20, ts21, ts22],
    ids=[20, 21, 22],
    symbols=["chr20", "chr21", "chr22"],
    types=["A", "A", "A"],  # Autosomes
    shared_nodes="samples",
)

# Compute statistics across all chromosomes
diversity = ta.stats.diversity()

# Save and load
ta.dump("my_genome_trees.zip")
ta_loaded = tmc.load("my_genome_trees.zip")
```

## Key Features

### TreesAssemblage

The central object holding a collection of per-contig tree sequences:

```python
# Access by symbol or id
ts = ta.contig("chr20")        # by symbol
ts = ta.contig(20)              # by id

# Iterate over contigs
for key, ts in ta.items():
    print(f"{key.symbol}: {ts.num_sites} sites")

# Properties
print(ta.total_sequence_length)     # Sum of all contig lengths
print(ta.num_contigs)               # Number of contigs
print(ta.global_phased_node_ids)    # Shared sample node IDs
```

### Subsetting and Reindexing

```python
# Subset by type, symbol, id, or index
autosomes = ta.subset(type="A")  # includes all of them in this case
chr20 = ta.subset(symbols=["chr20", "chr21"])

# Reindex contigs to 0..N
ta_reindexed = ta.reindex()

# Reorder contigs
ta_ordered = ta.reindex(order=["chr22", "chr21", "chr20"])
```

### Statistics

```python
# Across all contigs (requires samples to be globally phased)
diversity = ta.stats.diversity()

# With specific samples
diversity = ta.stats.diversity(sample_sets=[[0, 1, 2], [3, 4, 5]])

# By individual across all contigs
diversity = ta.stats.diversity(individuals=[0, 1])
```

### Format Conversion

```python
# From/to single tree sequence (merges all contigs)
ts_merged = ta.to_ts()
ta_from_merged = tmc.from_ts(ts_merged)

# From SLiM archives
ta_slim = tmc.from_slim(tree_sequences)
```

## Cross-Chromosome Phasing

The library tracks which sample nodes are **shared (phased) across all chromosomes** via the `IS_SHARED` flag:

- **Globally phased nodes**: Sample IDs that appear in all chromosomes with identical records
- **ARGs with nonglobal samples**: Archives where not all samples are phased across all chromosomes

Cross-chromosome stats automatically validate that sample IDs are globally phased to ensure correctness.

## Installation

```bash
pip install tskit_multichrom
```

Requires Python 3.9+, [tskit](https://tskit.dev/tskit), and [tszip](https://tskit.dev/tszip)

## Documentation

See [PROPOSAL.md](PROPOSAL.md) for detailed design specifications and implementation status.
