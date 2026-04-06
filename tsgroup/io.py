"""
Input/output functions for TreeSequenceGroup.

A trees archive can be stored in two formats:

1. **Directory** — a directory conventionally ending in ``_trees``, containing
   one ``.trees`` file (or ``.tsz`` if compressed with tszip) per contig.  The
   filenames are ``<symbol>.trees`` / ``<symbol>.tsz``.

2. **Zip file** — conventionally ending in ``_trees.zip``, an uncompressed
   (``ZIP_STORED``) zip of the same ``.trees`` files.  Compression is *not*
   applied inside the zip; use the tszip-based ``compress`` / ``decompress``
   helpers if you want to reduce file size.

In both cases each ``.trees`` file must carry a top-level ``'contig'`` metadata
key (see :class:`tskit.MetadataSchema.permissive_json`).
"""

import os
import shutil
import tempfile
import zipfile

import tskit
import tszip

from .core import ContigKey, TreeSequenceGroup, make_contig_key
from .flags import CONTIG_METADATA_KEY


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _load_tree_sequences_from_dir(path):
    """Load all .trees (or .tsz) files from a directory and return {key: ts}."""
    tree_sequences = {}
    entries = sorted(os.listdir(path))
    found = [
        e for e in entries if e.endswith(".trees") or e.endswith(".tsz")
    ]
    if not found:
        raise ValueError(f"No .trees or .tsz files found in directory {path!r}")

    for name in found:
        fpath = os.path.join(path, name)
        if name.endswith(".tsz"):
            ts = tszip.load(fpath)
        else:
            ts = tskit.load(fpath)

        meta = ts.metadata
        if not isinstance(meta, dict) or CONTIG_METADATA_KEY not in meta:
            raise ValueError(
                f"Tree sequence in {name!r} is missing "
                f"'{CONTIG_METADATA_KEY}' key in top-level metadata. "
                f"Got metadata: {meta!r}"
            )
        key = make_contig_key(meta[CONTIG_METADATA_KEY])
        if key in tree_sequences:
            raise ValueError(f"Duplicate ContigKey {key} found in archive")
        tree_sequences[key] = ts

    return tree_sequences


def _load_tree_sequences_from_zip(path):
    """Load all .trees (or .tsz) files from a zip archive and return {key: ts}."""
    tree_sequences = {}
    with zipfile.ZipFile(path, "r") as zf:
        names = [
            n for n in zf.namelist()
            if n.endswith(".trees") or n.endswith(".tsz")
        ]
        if not names:
            raise ValueError(f"No .trees or .tsz files found in {path!r}")

        with tempfile.TemporaryDirectory() as tmpdir:
            for name in names:
                dest = os.path.join(tmpdir, os.path.basename(name))
                with zf.open(name) as src, open(dest, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                if name.endswith(".tsz"):
                    ts = tszip.load(dest)
                else:
                    ts = tskit.load(dest)

                meta = ts.metadata
                if not isinstance(meta, dict) or CONTIG_METADATA_KEY not in meta:
                    raise ValueError(
                        f"Tree sequence in {name!r} is missing "
                        f"'{CONTIG_METADATA_KEY}' key in top-level metadata. "
                        f"Got metadata: {meta!r}"
                    )
                key = make_contig_key(meta[CONTIG_METADATA_KEY])
                if key in tree_sequences:
                    raise ValueError(
                        f"Duplicate ContigKey {key} found in archive"
                    )
                tree_sequences[key] = ts

    return tree_sequences


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load(path):
    """
    Load a :class:`~tsgroup.TreeSequenceGroup` from a trees archive.

    The archive may be:

    * A **directory** ending in ``_trees``, containing ``.trees`` or ``.tsz``
      files named ``<symbol>.trees`` / ``<symbol>.tsz``.
    * A **zip file** ending in ``_trees.zip``, containing ``.trees`` files.

    Parameters
    ----------
    path : str or path-like

    Returns
    -------
    TreeSequenceGroup
    """
    path = str(path)

    if os.path.isdir(path):
        tree_sequences = _load_tree_sequences_from_dir(path)
    elif zipfile.is_zipfile(path):
        tree_sequences = _load_tree_sequences_from_zip(path)
    else:
        raise ValueError(
            f"{path!r} is not a valid trees archive "
            "(expected a directory or a zip file)"
        )

    return TreeSequenceGroup(tree_sequences)


def dump(assemblage, path, *, compress=False):
    """
    Save a :class:`~tsgroup.TreeSequenceGroup` to a trees archive.

    Two output modes are available, selected by the *path* argument:

    * If *path* ends in ``_trees.zip`` (or ``.zip``), an **uncompressed zip**
      (``ZIP_STORED``) is written containing one ``<symbol>.trees`` file per
      contig.  Pass ``compress=True`` to use tszip inside the zip instead
      (``.tsz`` entries).
    * Otherwise *path* is treated as a **directory** and is created if it does
      not exist.  Each contig is written as ``<symbol>.trees`` (or
      ``<symbol>.tsz`` when ``compress=True``).

    Parameters
    ----------
    assemblage : TreeSequenceGroup
    path : str or path-like
        Destination path (directory or zip file).
    compress : bool
        If True, compress each tree sequence with tszip (produces ``.tsz`` files).
        Requires the ``tszip`` package.
    """
    path = str(path)
    path_lower = path.lower()
    is_zip = path_lower.endswith(".zip") or path_lower.endswith("_trees.zip")

    if is_zip:
        _dump_zip(assemblage, path, compress=compress)
    else:
        _dump_dir(assemblage, path, compress=compress)


def _dump_dir(assemblage, path, *, compress=False):
    """Write a trees archive to a directory."""
    os.makedirs(path, exist_ok=True)
    for key in assemblage.contigs:
        ts = assemblage[key]
        if compress:
            dest = os.path.join(path, f"{key.symbol}.tsz")
            tszip.compress(ts, dest)
        else:
            dest = os.path.join(path, f"{key.symbol}.trees")
            ts.dump(dest)


def _dump_zip(assemblage, path, *, compress=False):
    """Write a trees archive to an uncompressed zip file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as zf:
            for key in assemblage.contigs:
                ts = assemblage[key]
                if compress:
                    tmp_path = os.path.join(tmpdir, f"{key.symbol}.tsz")
                    tszip.compress(ts, tmp_path)
                    zf.write(tmp_path, arcname=f"{key.symbol}.tsz")
                else:
                    tmp_path = os.path.join(tmpdir, f"{key.symbol}.trees")
                    ts.dump(tmp_path)
                    zf.write(tmp_path, arcname=f"{key.symbol}.trees")
