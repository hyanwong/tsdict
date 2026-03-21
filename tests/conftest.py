"""
Test fixtures and helpers for tskit_multichrom tests.
"""

import numpy as np
import tskit
import tskit_multichrom as tmc


CONTIG_SCHEMA = {
    "codec": "json",
    "type": "object",
    "properties": {
        "contig": {
            "type": "object",
            "properties": {
                "index": {"type": "integer"},
                "id": {"type": "integer"},
                "symbol": {"type": "string"},
                "type": {"type": "string"},
            },
            "required": ["index", "id", "symbol", "type"],
        }
    },
    "additionalProperties": True,
}

NODE_SCHEMA = {
    "codec": "json",
    "type": "object",
    "properties": {
        "is_vacant": {"type": "integer"},
    },
    "additionalProperties": True,
}


def make_ts(
    seq_len=1000,
    num_samples=4,
    contig_meta=None,
    mark_shared=True,
):
    """
    Create a minimal tree sequence for testing.

    Parameters
    ----------
    seq_len : float
    num_samples : int  (must be even — pairs of haplotypes)
    contig_meta : dict with keys index, id, symbol, type (or None for no metadata)
    mark_shared : bool  — if True, set NODE_IS_SHARED on all nodes
    """
    assert num_samples % 2 == 0, "num_samples must be even"
    num_individuals = num_samples // 2

    tables = tskit.TableCollection(sequence_length=seq_len)
    tables.nodes.metadata_schema = tskit.MetadataSchema(NODE_SCHEMA)
    tables.populations.add_row()

    for _ in range(num_individuals):
        tables.individuals.add_row()

    node_flags = tskit.NODE_IS_SAMPLE
    if mark_shared:
        node_flags |= tmc.NODE_IS_SHARED

    for ind_id in range(num_individuals):
        tables.nodes.add_row(
            flags=node_flags,
            time=0,
            population=0,
            individual=ind_id,
            metadata={"is_vacant": 0},
        )
        tables.nodes.add_row(
            flags=node_flags,
            time=0,
            population=0,
            individual=ind_id,
            metadata={"is_vacant": 0},
        )

    # One ancestral node (not a sample, but shared too)
    anc_flags = tmc.NODE_IS_SHARED if mark_shared else 0
    tables.nodes.add_row(
        flags=anc_flags,
        time=1.0,
        population=0,
        metadata={"is_vacant": 0},
    )
    anc_id = num_samples

    for sample_id in range(num_samples):
        tables.edges.add_row(left=0, right=seq_len, parent=anc_id, child=sample_id)

    tables.sort()

    if contig_meta is not None:
        tables.metadata_schema = tskit.MetadataSchema(CONTIG_SCHEMA)
        tables.metadata = {"contig": contig_meta}

    return tables.tree_sequence()


def make_two_contig_archive(mark_shared=True, num_samples=4):
    """Return a TreesAssemblage with two autosomes."""
    ts1 = make_ts(
        seq_len=1000,
        num_samples=num_samples,
        contig_meta={"index": 0, "id": 0, "symbol": "chr1", "type": "A"},
        mark_shared=mark_shared,
    )
    ts2 = make_ts(
        seq_len=2000,
        num_samples=num_samples,
        contig_meta={"index": 1, "id": 1, "symbol": "chr2", "type": "A"},
        mark_shared=mark_shared,
    )
    return tmc.TreesAssemblage(
        {
            tmc.ContigKey(0, 0, "chr1", "A"): ts1,
            tmc.ContigKey(1, 1, "chr2", "A"): ts2,
        }
    )


def make_autosomes_plus_x_archive(num_samples=4):
    """
    Return an A/A/X TreesAssemblage where chrX has fewer sample nodes.

    chr1/chr2 have shared sample nodes 0..num_samples-1.
    chrX keeps only the first half as samples/shared; the rest are nonsample,
    nonshared. This gives different ``.samples()`` across chromosome types and
    exercises vacancy handling on roundtrip conversion.
    """
    ts1 = make_ts(
        seq_len=1000,
        num_samples=num_samples,
        contig_meta={"index": 0, "id": 0, "symbol": "chr1", "type": "A"},
        mark_shared=False,
    )
    ts2 = make_ts(
        seq_len=1200,
        num_samples=num_samples,
        contig_meta={"index": 1, "id": 1, "symbol": "chr2", "type": "A"},
        mark_shared=False,
    )
    tsx = make_ts(
        seq_len=800,
        num_samples=num_samples,
        contig_meta={"index": 2, "id": 2, "symbol": "chrX", "type": "X"},
        mark_shared=False,
    )

    result = {}
    for key, ts, mark_samples_shared in [
        (tmc.ContigKey(0, 0, "chr1", "A"), ts1, True),
        (tmc.ContigKey(1, 1, "chr2", "A"), ts2, True),
        (tmc.ContigKey(2, 2, "chrX", "X"), tsx, False),
    ]:
        tables = ts.dump_tables()
        if mark_samples_shared:
            flags = tables.nodes.flags.copy()
            sample_ids = ts.samples()
            flags[sample_ids] = flags[sample_ids] | tmc.NODE_IS_SHARED
            tables.nodes.flags = flags
        else:
            # chrX: keep only first half of samples as sample/shared
            flags = tables.nodes.flags.copy()
            sample_ids = ts.samples()
            keep_n = len(sample_ids) // 2
            keep_sample_ids = sample_ids[:keep_n]
            drop_sample_ids = sample_ids[keep_n:]

            # Keep sample status on first half and mark them shared.
            flags[keep_sample_ids] = flags[keep_sample_ids] | tmc.NODE_IS_SHARED

            # Remove sample + shared flags from second half.
            clear_bits = np.uint32(tskit.NODE_IS_SAMPLE | tmc.NODE_IS_SHARED)
            flags[drop_sample_ids] = flags[drop_sample_ids] & (~clear_bits)
            tables.nodes.flags = flags
        result[key] = tables.tree_sequence()

    return tmc.TreesAssemblage(result)
