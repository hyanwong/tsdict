"""
tsgroup: A wrapper on top of tskit for storing and analysing
multiple chromosomes (contigs).

The central object is :class:`TreeSequenceGroup` (abbreviated ``tsg``), which holds
a collection of :class:`tskit.TreeSequence` objects — one per contig — stored
in a dictionary keyed by :class:`ContigKey` namedtuples.

Quick start::

    import tsgroup

    # Load from a directory or zip archive
    tsg = tsgroup.load("genome_trees")

    # Create from a list of tree sequences
    tsg = tsgroup.from_tree_sequences(
        [ts1, ts2],
        ids=[1, 2],
        symbols=["chr1", "chr2"],
        types=["A", "A"],
    )

    # Access a contig by symbol
    chr1_ts = tsg.contig("chr1")

    # Save back to a directory
    tsg.dump("genome_trees")

    # Convert to a single tree sequence
    ts = tsg.to_ts()

    # Convert back
    tsg2 = tsgroup.from_ts(ts)
"""

from .core import (
    ContigKey,
    TreeSequenceGroup,
    make_contig_key,
)
from .convert import from_slim, from_tree_sequences, from_ts
from .flags import ARCHIVE_EXTENSION, CONTIG_METADATA_KEY, NODE_IS_SHARED
from .io import load
from ._version import tsgroup_version as __version__

__all__ = [
    # Core
    "TreeSequenceGroup",
    "ContigKey",
    "make_contig_key",
    # I/O
    "load",
    # Conversion
    "from_ts",
    "from_slim",
    "from_tree_sequences",
    # Flags
    "NODE_IS_SHARED",
    "CONTIG_METADATA_KEY",
    "ARCHIVE_EXTENSION",
]

