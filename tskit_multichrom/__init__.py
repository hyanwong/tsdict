"""
tskit_multichrom: A wrapper on top of tskit for storing and analysing
multiple chromosomes (contigs).

The central object is :class:`TreesAssemblage` (abbreviated ``ta``), which holds
a collection of :class:`tskit.TreeSequence` objects — one per contig — stored
in a dictionary keyed by :class:`ContigKey` namedtuples.

Quick start::

    import tskit_multichrom as tmc

    # Load from a directory or zip archive
    ta = tmc.load("genome_trees")

    # Create from a list of tree sequences
    ta = tmc.from_tree_sequences(
        [ts1, ts2],
        ids=[1, 2],
        symbols=["chr1", "chr2"],
        types=["A", "A"],
    )

    # Access a contig by symbol
    chr1_ts = ta.contig("chr1")

    # Save back to a directory
    ta.dump("genome_trees")
    # or equivalently
    tmc.dump(ta, "genome_trees")

    # Convert to a single tree sequence
    ts = tmc.to_ts(ta)

    # Convert back
    ta2 = tmc.from_ts(ts)
"""

from .core import (
    ContigKey,
    TreesAssemblage,
    make_contig_key,
)
from .convert import from_slim, from_tree_sequences, from_ts, to_ts
from .flags import ARCHIVE_EXTENSION, CONTIG_METADATA_KEY, NODE_IS_SHARED
from .io import dump, load
from ._version import tskit_multichrom_version as __version__

__all__ = [
    # Core
    "TreesAssemblage",
    "ContigKey",
    "make_contig_key",
    # I/O
    "load",
    "dump",
    # Conversion
    "to_ts",
    "from_ts",
    "from_slim",
    "from_tree_sequences",
    # Flags
    "NODE_IS_SHARED",
    "CONTIG_METADATA_KEY",
    "ARCHIVE_EXTENSION",
]

