#!/usr/bin/env python3
import tskit_multichrom as tmc
from tests.conftest import make_ts

# Create a simple tree sequence without contig metadata
ts = make_ts(seq_len=1000, num_samples=4, mark_shared=False)

# 1. Default behavior: no nodes marked as IS_SHARED
print("1. Default (shared_nodes=None):")
ta1 = tmc.from_tree_sequences(
    [ts],
    ids=[1],
    symbols=["chr1"],
    types=["A"],
)
ts_result = ta1.contig("chr1")
is_shared_count = sum(1 for nid in range(ts_result.num_nodes) 
                      if ts_result.tables.nodes[nid].flags & tmc.NODE_IS_SHARED)
print(f"   Nodes marked IS_SHARED: {is_shared_count} (expected 0)")
assert is_shared_count == 0, "Default should not mark any nodes"

# 2. Explicit shared_nodes="samples"
print("\n2. With shared_nodes='samples':")
ta2 = tmc.from_tree_sequences(
    [ts],
    ids=[1],
    symbols=["chr1"],
    types=["A"],
    shared_nodes="samples",
)
ts_result = ta2.contig("chr1")
is_shared_count = sum(1 for nid in range(ts_result.num_nodes) 
                      if ts_result.tables.nodes[nid].flags & tmc.NODE_IS_SHARED)
sample_count = len(ts_result.samples())
print(f"   Nodes marked IS_SHARED: {is_shared_count} (expected {sample_count})")
assert is_shared_count == sample_count, f"Expected {sample_count} nodes marked"

# 3. Explicit list of node IDs
print("\n3. With shared_nodes=[0, 1]:")
ta3 = tmc.from_tree_sequences(
    [ts],
    ids=[1],
    symbols=["chr1"],
    types=["A"],
    shared_nodes=[0, 1],
)
ts_result = ta3.contig("chr1")
is_shared_count = sum(1 for nid in range(ts_result.num_nodes) 
                      if ts_result.tables.nodes[nid].flags & tmc.NODE_IS_SHARED)
print(f"   Nodes marked IS_SHARED: {is_shared_count} (expected 2)")
assert is_shared_count == 2, "Expected 2 nodes marked"
assert ts_result.tables.nodes[0].flags & tmc.NODE_IS_SHARED
assert ts_result.tables.nodes[1].flags & tmc.NODE_IS_SHARED

print("\n✓ All three usage patterns work correctly!")
