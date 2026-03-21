"""
Tests for TreesArchive I/O (load/dump).
"""

import os
import tempfile

import pytest
import tskit

import tskit_multichrom as tmc
from tests.conftest import make_two_contig_archive


class TestDumpLoad:
    def test_roundtrip(self):
        ta = make_two_contig_archive()
        with tempfile.NamedTemporaryFile(suffix=".tsa", delete=False) as f:
            path = f.name
        try:
            tmc.dump(ta, path)
            ta2 = tmc.load(path)
            assert ta2.num_contigs == 2
            assert ta2.total_sequence_length == ta.total_sequence_length
        finally:
            os.unlink(path)

    def test_roundtrip_contigs(self):
        ta = make_two_contig_archive()
        with tempfile.NamedTemporaryFile(suffix=".tsa", delete=False) as f:
            path = f.name
        try:
            tmc.dump(ta, path)
            ta2 = tmc.load(path)
            assert ta2.chr("chr1").sequence_length == 1000
            assert ta2.chr("chr2").sequence_length == 2000
        finally:
            os.unlink(path)

    def test_load_non_zip_raises(self):
        with tempfile.NamedTemporaryFile(suffix=".tsa", delete=False) as f:
            f.write(b"not a zip file")
            path = f.name
        try:
            with pytest.raises(ValueError, match="valid trees archive"):
                tmc.load(path)
        finally:
            os.unlink(path)

    def test_load_empty_zip_raises(self):
        import zipfile

        with tempfile.NamedTemporaryFile(suffix=".tsa", delete=False) as f:
            path = f.name
        try:
            with zipfile.ZipFile(path, "w"):
                pass
            with pytest.raises(ValueError, match="No .trees files"):
                tmc.load(path)
        finally:
            os.unlink(path)

    def test_node_flags_preserved(self):
        ta = make_two_contig_archive(mark_shared=True)
        with tempfile.NamedTemporaryFile(suffix=".tsa", delete=False) as f:
            path = f.name
        try:
            tmc.dump(ta, path)
            ta2 = tmc.load(path)
            ts = ta2.chr("chr1")
            # All nodes should have IS_SHARED flag set
            for node_id in range(ts.num_nodes):
                flags = ts.tables.nodes[node_id].flags
                assert flags & tmc.NODE_IS_SHARED, (
                    f"Node {node_id} should have IS_SHARED set, got flags={flags}"
                )
        finally:
            os.unlink(path)
