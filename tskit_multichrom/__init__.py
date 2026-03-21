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
    make_permissive_contig_schema,
)
from .convert import from_slim, from_ts, to_ts
from .flags import ARCHIVE_EXTENSION, CONTIG_METADATA_KEY, NODE_IS_SHARED
from .io import dump, load

__all__ = [
    # Core
    "TreesAssemblage",
    "ContigKey",
    "make_contig_key",
    "make_permissive_contig_schema",
    # I/O
    "load",
    "dump",
    # Conversion
    "to_ts",
    "from_ts",
    "from_slim",
    # Flags
    "NODE_IS_SHARED",
    "CONTIG_METADATA_KEY",
    "ARCHIVE_EXTENSION",
]

__version__ = "0.1.0"
