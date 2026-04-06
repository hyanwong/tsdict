# Usage

## Quick Start

```python
import msprime
import tsgroup

# Simulate two chromosomes with msprime
ts20 = msprime.simulate(sample_size=8, length=1e6, recombination_rate=1e-8, random_seed=1)
ts21 = msprime.simulate(sample_size=8, length=5e5, recombination_rate=1e-8, random_seed=2)

# Assemble into a multi-chromosome object
# shared_nodes="samples" marks all sample nodes as cross-chromosome phased
tsg = tsgroup.from_tree_sequences(
    [ts20, ts21],
    ids=[20, 21],
    symbols=["chr20", "chr21"],
    types=["A", "A"],
    shared_nodes="samples",
)

print(tsg)
# TreeSequenceGroup(contigs=['chr20', 'chr21'], total_length=1500000.0)
```

## Saving and Loading

Archives can be stored as a **directory** (conventionally `*_trees`) or as a
**zip file** (conventionally `*_trees.zip`):

```python
# Save to zip
tsg.dump("my_genome_trees.zip")

# Load back
tsg2 = tsgroup.load("my_genome_trees.zip")

# Compress individual tree sequences with tszip
tsg.dump("my_genome_trees.zip", compress=True)
```

## Accessing Contigs

```python
# By symbol
ts = tsg.contig("chr20")

# By integer id
ts = tsg.contig(20)

# Iterate over all contigs
for key, ts in tsg.items():
    print(f"{key.symbol}: {ts.num_sites} sites, length {ts.sequence_length}")
```

## Key Properties

```python
tsg.num_contigs                # number of contigs
tsg.total_sequence_length      # sum of all contig lengths
tsg.global_phased_node_ids     # sample nodes shared across all contigs
tsg.shared_node_ids            # nodes with IS_SHARED in at least one contig
tsg.is_nonglobal_sample_arg    # True if some samples are not globally phased
```

## Subsetting

Subsetting is cheap — it does not copy tree sequences:

```python
# Keep only autosomes
autosomes = tsg.subset(type="A")

# Keep specific chromosomes by symbol
subset = tsg.subset(symbols=["chr20", "chr21"])

# Keep by integer contig id
subset = tsg.subset(ids=[20, 22])
```

## Cross-chromosome Statistics

Statistics require that the relevant sample nodes are globally phased (present
and marked as {attr}`~tsgroup.NODE_IS_SHARED` in **all** contigs).

```python
# Diversity across all contigs (weighted by span)
pi = tsg.stats.diversity()

# Diversity for specific sample sets
pi = tsg.stats.diversity(sample_sets=[[0, 2, 4], [1, 3, 5]])

# PCA across all contigs
result = tsg.stats.pca(num_components=5)
print(result.coordinates)  # shape (n_samples, 5)
```

## Converting to/from a Single Tree Sequence

```python
# Merge all contigs into one tree sequence (coordinates concatenated)
ts = tsg.to_ts()

# Split back into a TreeSequenceGroup
tsg2 = tsgroup.from_ts(ts)
```

## SLiM Compatibility

SLiM (prior to v5) produces one tree sequence per chromosome where every node
is implicitly shared.  Use {func}`~tsgroup.from_slim` to import them:

```python
# slim_ts_list is a list of SLiM-generated tskit.TreeSequence objects
tsg = tsgroup.from_slim(slim_ts_list)
```

## Simplifying

```python
# Simplify to globally phased sample nodes only
simplified = tsg.simplify()

# Simplify to a specific set of sample nodes
simplified = tsg.simplify(samples=[0, 2, 4])

# Simplify per-individual (works even with non-globally phased nodes)
simplified = tsg.simplify(individuals=[0, 1, 2])
```
