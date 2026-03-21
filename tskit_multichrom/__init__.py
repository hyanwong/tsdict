"""
tskit_multichrom: A wrapper on top of tskit for storing and analysing
multiple chromosomes (contigs).

The central object is :class:`TreesArchive` (abbreviated ``ta``), which holds
a collection of :class:`tskit.TreeSequence` objects — one per contig — stored
in a dictionary keyed by :class:`ContigKey` namedtuples.

Quick start::

    import tskit_multichrom as tmc

    # Load from a .tsa archive file
    ta = tmc.load("genome.tsa")

    # Access a contig by symbol
    chr1_ts = ta.chr("chr1")

    # Save back to file
    tmc.dump(ta, "genome.tsa")

    # Convert to a single tree sequence
    ts = tmc.to_tree_sequence(ta)

    # Convert back
    ta2 = tmc.from_tree_sequence(ts)
"""

from .core import ContigKey, TreesArchive, make_contig_key, make_permissive_contig_schema
from .convert import from_slim, from_tree_sequence, to_tree_sequence
from .flags import ARCHIVE_EXTENSION, CONTIG_METADATA_KEY, NODE_IS_SHARED
from .io import dump, load

__all__ = [
    # Core
    "TreesArchive",
    "ContigKey",
    "make_contig_key",
    "make_permissive_contig_schema",
    # I/O
    "load",
    "dump",
    # Conversion
    "to_tree_sequence",
    "from_tree_sequence",
    "from_slim",
    # Flags
    "NODE_IS_SHARED",
    "CONTIG_METADATA_KEY",
    "ARCHIVE_EXTENSION",
]

__version__ = "0.1.0"
