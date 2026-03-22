"""tsdict public API."""

from .core import ContigKey, TreeSequenceDictionary, make_contig_key
from .convert import from_slim, from_tree_sequences, from_ts
from .flags import ARCHIVE_EXTENSION, CONTIG_METADATA_KEY, NODE_IS_SHARED
from .io import load

from ._version import tsdict_version as __version__

TreeSequenceDictionary = TreeSequenceDictionary
TSDict = TreeSequenceDictionary

__all__ = [
    "TreeSequenceDictionary",
    "TSDict",
    "TreeSequenceDictionary",
    "ContigKey",
    "make_contig_key",
    "load",
    "from_ts",
    "from_slim",
    "from_tree_sequences",
    "NODE_IS_SHARED",
    "CONTIG_METADATA_KEY",
    "ARCHIVE_EXTENSION",
]
