"""
Tests for TreeSequenceGroup I/O (load/dump).
"""

import os
import tempfile
import zipfile

import pytest
import tskit

import tsgroup
from tests.conftest import make_two_contig_archive


class TestDumpLoadZip:
    """Tests for the _trees.zip format."""

    def test_roundtrip(self):
        tsg = make_two_contig_archive()
        with tempfile.NamedTemporaryFile(suffix="_trees.zip", delete=False) as f:
            path = f.name
        try:
            tsg.dump(path)
            tsg2 = tsgroup.load(path)
            assert tsg2.num_contigs == 2
            assert tsg2.total_sequence_length == tsg.total_sequence_length
        finally:
            os.unlink(path)

    def test_roundtrip_contigs(self):
        tsg = make_two_contig_archive()
        with tempfile.NamedTemporaryFile(suffix="_trees.zip", delete=False) as f:
            path = f.name
        try:
            tsg.dump(path)
            tsg2 = tsgroup.load(path)
            assert tsg2.contig("chr1").sequence_length == 1000
            assert tsg2.contig("chr2").sequence_length == 2000
        finally:
            os.unlink(path)

    def test_filenames_use_symbol(self):
        """Files inside the zip should use symbol, not index."""
        tsg = make_two_contig_archive()
        with tempfile.NamedTemporaryFile(suffix="_trees.zip", delete=False) as f:
            path = f.name
        try:
            tsg.dump(path)
            with zipfile.ZipFile(path, "r") as zf:
                names = zf.namelist()
            assert "chr1.trees" in names
            assert "chr2.trees" in names
        finally:
            os.unlink(path)

    def test_zip_is_stored_not_compressed(self):
        """The zip should use ZIP_STORED (no deflate)."""
        tsg = make_two_contig_archive()
        with tempfile.NamedTemporaryFile(suffix="_trees.zip", delete=False) as f:
            path = f.name
        try:
            tsg.dump(path)
            with zipfile.ZipFile(path, "r") as zf:
                for info in zf.infolist():
                    assert info.compress_type == zipfile.ZIP_STORED, (
                        f"{info.filename} should be stored, not compressed"
                    )
        finally:
            os.unlink(path)

    def test_load_non_zip_raises(self):
        with tempfile.NamedTemporaryFile(suffix="_trees.zip", delete=False) as f:
            f.write(b"not a zip file")
            path = f.name
        try:
            with pytest.raises(ValueError, match="valid trees archive"):
                tsgroup.load(path)
        finally:
            os.unlink(path)

    def test_load_empty_zip_raises(self):
        with tempfile.NamedTemporaryFile(suffix="_trees.zip", delete=False) as f:
            path = f.name
        try:
            with zipfile.ZipFile(path, "w"):
                pass
            with pytest.raises(ValueError, match="No .trees or .tsz files"):
                tsgroup.load(path)
        finally:
            os.unlink(path)

    def test_node_flags_preserved(self):
        tsg = make_two_contig_archive(mark_shared=True)
        with tempfile.NamedTemporaryFile(suffix="_trees.zip", delete=False) as f:
            path = f.name
        try:
            tsg.dump(path)
            tsg2 = tsgroup.load(path)
            ts = tsg2.contig("chr1")
            for node_id in range(ts.num_nodes):
                flags = ts.tables.nodes[node_id].flags
                assert flags & tsgroup.NODE_IS_SHARED, (
                    f"Node {node_id} should have IS_SHARED set, got flags={flags}"
                )
        finally:
            os.unlink(path)


class TestDumpLoadDir:
    """Tests for the directory (_trees) format."""

    def test_roundtrip(self):
        tsg = make_two_contig_archive()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "genome_trees")
            tsg.dump(path)
            tsg2 = tsgroup.load(path)
            assert tsg2.num_contigs == 2
            assert tsg2.total_sequence_length == tsg.total_sequence_length

    def test_roundtrip_contigs(self):
        tsg = make_two_contig_archive()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "genome_trees")
            tsg.dump(path)
            tsg2 = tsgroup.load(path)
            assert tsg2.contig("chr1").sequence_length == 1000
            assert tsg2.contig("chr2").sequence_length == 2000

    def test_filenames_use_symbol(self):
        tsg = make_two_contig_archive()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "genome_trees")
            tsg.dump(path)
            files = os.listdir(path)
        assert "chr1.trees" in files
        assert "chr2.trees" in files

    def test_load_empty_dir_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(ValueError, match="No .trees"):
                tsgroup.load(tmpdir)


class TestDumpMethod:
    """Tests for TreeSequenceGroup.dump() instance method."""

    def test_dump_method_zip(self):
        tsg = make_two_contig_archive()
        with tempfile.NamedTemporaryFile(suffix="_trees.zip", delete=False) as f:
            path = f.name
        try:
            tsg.dump(path)
            tsg2 = tsgroup.load(path)
            assert tsg2.num_contigs == 2
            assert tsg2.total_sequence_length == tsg.total_sequence_length
        finally:
            os.unlink(path)

    def test_dump_method_dir(self):
        tsg = make_two_contig_archive()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "genome_trees")
            tsg.dump(path)
            tsg2 = tsgroup.load(path)
            assert tsg2.num_contigs == 2
