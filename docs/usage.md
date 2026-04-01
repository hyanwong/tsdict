# Usage

## Quick Start

```python
import msprime
import tskit_multichrom as tmc

# Simulate two chromosomes with msprime
ts20 = msprime.simulate(sample_size=8, length=1e6, recombination_rate=1e-8, random_seed=1)
ts21 = msprime.simulate(sample_size=8, length=5e5, recombination_rate=1e-8, random_seed=2)

# Assemble into a multi-chromosome object
# shared_nodes="samples" marks all sample nodes as cross-chromosome phased
tsd = tmc.from_tree_sequences(
    [ts20, ts21],
    ids=[20, 21],
    symbols=["chr20", "chr21"],
    types=["A", "A"],
    shared_nodes="samples",
)

print(tsd)
# TreeSequenceDictionary(contigs=['chr20', 'chr21'], total_length=1500000.0)
```

## Saving and Loading

Archives can be stored as a **directory** (conventionally `*_trees`) or as a
**zip file** (conventionally `*_trees.zip`):

```python
# Save to zip
tsd.dump("my_genome_trees.zip")

# Load back
tsd2 = tmc.load("my_genome_trees.zip")

# Compress individual tree sequences with tszip
tsd.dump("my_genome_trees.zip", compress=True)
```

## Accessing Contigs

```python
# By symbol
ts = tsd.contig("chr20")

# By integer id
ts = tsd.contig(20)

# Iterate over all contigs
for key, ts in tsd.items():
    print(f"{key.symbol}: {ts.num_sites} sites, length {ts.sequence_length}")
```

## Key Properties

```python
tsd.num_contigs                # number of contigs
tsd.total_sequence_length      # sum of all contig lengths
tsd.global_phased_node_ids     # sample nodes shared across all contigs
tsd.shared_node_ids            # nodes with IS_SHARED in at least one contig
tsd.is_nonglobal_sample_arg    # True if some samples are not globally phased
```

## Subsetting

Subsetting is cheap — it does not copy tree sequences:

```python
# Keep only autosomes
autosomes = tsd.subset(type="A")

# Keep specific chromosomes by symbol
subset = tsd.subset(symbols=["chr20", "chr21"])

# Keep by integer contig id
subset = tsd.subset(ids=[20, 22])
```

## Cross-chromosome Statistics

Statistics require that the relevant sample nodes are globally phased (present
and marked as {attr}`~tskit_multichrom.NODE_IS_SHARED` in **all** contigs).

```python
# Diversity across all contigs (weighted by span)
pi = tsd.stats.diversity()

# Diversity for specific sample sets
pi = tsd.stats.diversity(sample_sets=[[0, 2, 4], [1, 3, 5]])

# PCA across all contigs
result = tsd.stats.pca(num_components=5)
print(result.coordinates)  # shape (n_samples, 5)
```

## Converting to/from a Single Tree Sequence

```python
# Merge all contigs into one tree sequence (coordinates concatenated)
ts = tsd.to_ts()

# Split back into a TreeSequenceDictionary
tsd2 = tmc.from_ts(ts)
```

## SLiM Compatibility

SLiM (prior to v5) produces one tree sequence per chromosome where every node
is implicitly shared.  Use {func}`~tskit_multichrom.from_slim` to import them:

```python
# slim_ts_list is a list of SLiM-generated tskit.TreeSequence objects
tsd = tmc.from_slim(slim_ts_list)
```

## Simplifying

```python
# Simplify to globally phased sample nodes only
simplified = tsd.simplify()

# Simplify to a specific set of sample nodes
simplified = tsd.simplify(samples=[0, 2, 4])

# Simplify per-individual (works even with non-globally phased nodes)
simplified = tsd.simplify(individuals=[0, 1, 2])
```
