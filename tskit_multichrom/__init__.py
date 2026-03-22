"""
tskit_multichrom: A wrapper on top of tskit for storing and analysing
multiple chromosomes (contigs).

The central object is :class:`TreeSequenceDictionary` (abbreviated ``tsd``), which holds
a collection of :class:`tskit.TreeSequence` objects — one per contig — stored
in a dictionary keyed by :class:`ContigKey` namedtuples.

Quick start::

    import tskit_multichrom as tmc

    # Load from a directory or zip archive
    tsd = tmc.load("genome_trees")

    # Create from a list of tree sequences
    tsd = tmc.from_tree_sequences(
        [ts1, ts2],
        ids=[1, 2],
        symbols=["chr1", "chr2"],
        types=["A", "A"],
    )

    # Access a contig by symbol
    chr1_ts = tsd.contig("chr1")

    # Save back to a directory
    tsd.dump("genome_trees")

    # Convert to a single tree sequence
    ts = tsd.to_ts()

    # Convert back
    tsd2 = tmc.from_ts(ts)
"""

from .core import (
    ContigKey,
    TreeSequenceDictionary,
    make_contig_key,
)
from .convert import from_slim, from_tree_sequences, from_ts
from .flags import ARCHIVE_EXTENSION, CONTIG_METADATA_KEY, NODE_IS_SHARED
from .io import load
from ._version import tskit_multichrom_version as __version__

__all__ = [
    # Core
    "TreeSequenceDictionary",
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

