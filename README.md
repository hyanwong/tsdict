# tsdict

A Python library for efficiently storing and analyzing multiple chromosomes (contigs) using [tskit](https://tskit.dev). It allows you to work with a collection of tree sequences — one per chromosome — as a unified dictionary object while maintaining cross-chromosome sample identity.

## Why tsdict?

Population genetics simulators like **msprime** and **SLiM** (prior to v5) typically produce separate tree sequences for each chromosome. `tsdict` provides:

- **Unified storage**: Manage multiple per-chromosome tree sequences as a single `TreeSequenceDictionary` object
- **Cross-chromosome analysis**: Compute statistics across all chromosomes when samples are shared (e.g., diversity, PCA)
- **Efficient subsetting**: Create subsets (e.g., autosomes only) without copying underlying tree sequences
- **Format conversion**: Convert between multi-contig archives and single merged tree sequences
- **SLiM compatibility**: Seamless conversion to/from SLiM's tree sequence archives

## Quick Start

```python
import stdpopsim
import tsdict

# Simulate multiple chromosomes using stdpopsim
species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model("OutOfAfrica_2T12")
chr20 = species.get_contig("chr20", mutation_rate=model.mutation_rate)
chr21 = species.get_contig("chr21", mutation_rate=model.mutation_rate)
chr22 = species.get_contig("chr22", mutation_rate=model.mutation_rate)

samples = {'AFR': 3, 'EUR': 4}
engine = stdpopsim.get_engine("msprime")
ts20 = engine.simulate(model, chr20, samples)
ts21 = engine.simulate(model, chr21, samples)
ts22 = engine.simulate(model, chr22, samples)
first_8_haplosomes = ts20.samples()[0:8]

# Combine into a multi-chromosome assemblage
# shared_nodes can be a list of sample nodes to mark as shared across chromosomes,
# (required for cross-chromosome statistics) or "samples" to mark all samples as such
tsd = tsdict.from_tree_sequences(
    [ts20, ts21, ts22],
    ids=[20, 21, 22],
    symbols=["chr20", "chr21", "chr22"],
    types=["A", "A", "A"],  # Autosomes
    shared_nodes=first_8_haplosomes,  # Only mark some samples as across-chromosome phased
)

# Compute statistics across all chromosomes
tsd.stats.diversity()  # fails as not all samples are cross-phased
diversity = tsd.stats.diversity(first_8_haplosomes)  # passes
# Or simplify down to just the cross-phased samples
cross_phased_tsd = tsd.simplify(samples=first_8_haplosomes)
assert diversity == cross_phased_tsd.stats.diversity()

# Save and load
tsd.dump("my_genome_trees.zip")
tsd_loaded = tsdict.load("my_genome_trees.zip")
```

## Key Features

### TreeSequenceDictionary

The central object holding a collection of per-contig tree sequences:

```python
# Access by symbol or id
ts = tsd.contig("chr20")        # by symbol
ts = tsd.contig(20)             # by id

# Iterate over contigs
for key, ts in tsd.items():
    print(f"{key.symbol}: {ts.num_sites} sites")

# Properties
print(tsd.total_sequence_length)     # Sum of all contig lengths
print(tsd.num_contigs)               # Number of contigs
print(tsd.global_phased_node_ids)    # Shared sample node IDs
```

### Subsetting and Reindexing

```python
# Subset by type, symbol, id, or index
autosomes = tsd.subset(type="A")  # includes all the chromosomes in this case
chr21and22 = tsd.subset(symbols=["chr21", "chr22"])

# Reindex contigs to 0..N
tsd21_22_reindexed = chr21and22.reindex()

# Reorder contigs
tsd_reordered = tsd.reindex(order=["chr22", "chr21", "chr20"])
```

### Statistics

```python
# Across all contigs (requires samples to be globally phased)
diversity = cross_phased_tsd.stats.diversity()

# With specific samples
diversity_array = cross_phased_tsd.stats.diversity(sample_sets=[[0, 1, 2], [3, 4, 5]])

# Some stats can be run by individual (e.g. PCA)
# These do not require phasing across chromosomes
pca_result = tsd.stats.pca(num_components=2, individuals=[0, 1, 2, 3])
```

### Format Conversion

```python
# From/to single tree sequence (merges all contigs)
ts_merged = tsd.to_ts()
tsd_from_merged = tsdict.from_ts(ts_merged)

# From SLiM archives
tsd_slim = tsdict.from_slim(tree_sequences)
```

## Cross-Chromosome Phasing

The library tracks which sample nodes are **shared (phased) across all chromosomes** via the `IS_SHARED` flag:

- **Globally phased nodes**: Sample IDs that appear in all chromosomes with identical records
- **ARGs with nonglobal samples**: Archives where not all samples are phased across all chromosomes

Cross-chromosome stats automatically validate that sample IDs are globally phased to ensure correctness.

## Installation

For the moment, the package is in test mode and only available through git:

```bash
pip install git+https://github.com/hyanwong/tsdict
```

Requires Python 3.9+, [tskit](https://tskit.dev/tskit), and [tszip](https://tskit.dev/tszip)

## Documentation

See [PROPOSAL.md](PROPOSAL.md) for detailed design specifications and implementation status.
