"""
Flag constants for tsgroup.
"""

# Node flag indicating that this node is shared (cross-phased) across all
# contigs in a TreesArchive. This bit is set in the node flags field.
# We use bit 20 to avoid conflicts with tskit (bit 0) and SLiM (bits 16–19).
NODE_IS_SHARED = 1 << 20

# Default file suffix for an uncompressed trees archive zip file.
# A trees archive can also be a directory ending in ``_trees``.
ARCHIVE_EXTENSION = "_trees.zip"

# Name of the metadata key holding contig information in each tree sequence
CONTIG_METADATA_KEY = "contig"
