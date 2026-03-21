"""
Input/output functions for TreesArchive.

The trees archive file format is a zip file (extension ``.tsa``) containing
one ``.trees`` file per contig, named by index (``0.trees``, ``1.trees``, …).
"""

import os
import shutil
import tempfile
import zipfile

import tskit

from .core import ContigKey, TreesArchive, make_contig_key
from .flags import ARCHIVE_EXTENSION, CONTIG_METADATA_KEY


def load(path):
    """
    Load a :class:`TreesArchive` from a ``.tsa`` zip archive.

    Each constituent ``.trees`` file must have a top-level metadata key
    ``'contig'`` (see :func:`~tskit_multichrom.core.make_permissive_contig_schema`).

    Parameters
    ----------
    path : str or path-like
        Path to the ``.tsa`` file.

    Returns
    -------
    TreesArchive
    """
    path = str(path)
    if not zipfile.is_zipfile(path):
        raise ValueError(f"{path!r} is not a valid trees archive (zip) file")

    tree_sequences = {}

    with zipfile.ZipFile(path, "r") as zf:
        names = [n for n in zf.namelist() if n.endswith(".trees")]
        if not names:
            raise ValueError(f"No .trees files found in {path!r}")

        with tempfile.TemporaryDirectory() as tmpdir:
            for name in names:
                # Extract to a temp file and load
                dest = os.path.join(tmpdir, os.path.basename(name))
                with zf.open(name) as src, open(dest, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                ts = tskit.load(dest)

                # Extract contig key from metadata
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

    return TreesArchive(tree_sequences)


def dump(archive, path):
    """
    Save a :class:`TreesArchive` to a ``.tsa`` zip archive.

    Parameters
    ----------
    archive : TreesArchive
        The archive to save.
    path : str or path-like
        Destination path (should end in ``.tsa``).
    """
    path = str(path)

    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for key in archive.contigs:
                ts = archive[key]
                tmp_path = os.path.join(tmpdir, f"{key.index}.trees")
                ts.dump(tmp_path)
                zf.write(tmp_path, arcname=f"{key.index}.trees")
