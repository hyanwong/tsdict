"""
Tests for conversion between TreesAssemblage and a single TreeSequence.
"""

import pytest
import tskit

import tskit_multichrom as tmc
from tests.conftest import make_ts, make_two_contig_archive


class TestToTreeSequence:
    def test_basic(self):
        ta = make_two_contig_archive()
        ts = tmc.to_ts(ta)
        assert isinstance(ts, tskit.TreeSequence)
        assert ts.sequence_length == ta.total_sequence_length

    def test_metadata_present(self):
        ta = make_two_contig_archive()
        ts = tmc.to_ts(ta)
        assert "tskit_multichrom_contigs" in ts.metadata

    def test_contigs_in_metadata(self):
        ta = make_two_contig_archive()
        ts = tmc.to_ts(ta)
        meta = ts.metadata["tskit_multichrom_contigs"]
        assert len(meta) == 2
        symbols = [c["contig"]["symbol"] for c in meta]
        assert "chr1" in symbols
        assert "chr2" in symbols

    def test_empty_raises(self):
        ta = tmc.TreesAssemblage({})
        with pytest.raises(ValueError, match="empty"):
            tmc.to_ts(ta)

    def test_site_schema_mismatch_raises(self):
        ts1 = make_ts(
            seq_len=1000,
            contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"},
        )
        ts2 = make_ts(
            seq_len=2000,
            contig_meta={"index": 1, "id": 1, "symbol": "c2", "type": "A"},
        )
        # Give ts2 a different site schema
        tables = ts2.dump_tables()
        tables.sites.metadata_schema = tskit.MetadataSchema(
            {"codec": "json", "type": "object", "properties": {"x": {"type": "integer"}}}
        )
        ts2_mod = tables.tree_sequence()
        ta = tmc.TreesAssemblage(
            {
                tmc.ContigKey(0, 0, "c1", "A"): ts1,
                tmc.ContigKey(1, 1, "c2", "A"): ts2_mod,
            }
        )
        with pytest.raises(ValueError, match="[Ss]ite"):
            tmc.to_ts(ta)

    def test_shared_nodes_have_same_id(self):
        """Shared nodes should retain their IDs in the merged tree sequence."""
        ta = make_two_contig_archive(mark_shared=True)
        ts = tmc.to_ts(ta)

        # Nodes 0..4 (4 samples + 1 ancestral) should all have IS_SHARED set
        for node_id in range(5):
            flags = ts.tables.nodes[node_id].flags
            assert flags & tmc.NODE_IS_SHARED

    def test_sequence_length(self):
        ta = make_two_contig_archive()
        ts = tmc.to_ts(ta)
        assert ts.sequence_length == 3000  # 1000 + 2000

    def test_sites_are_shifted(self):
        """Sites from chr2 should be shifted by chr1's sequence_length."""
        ts1 = make_ts(
            seq_len=1000,
            contig_meta={"index": 0, "id": 0, "symbol": "chr1", "type": "A"},
        )
        ts2 = make_ts(
            seq_len=2000,
            contig_meta={"index": 1, "id": 1, "symbol": "chr2", "type": "A"},
        )
        # Add a site to each
        tables1 = ts1.dump_tables()
        tables1.sites.add_row(position=500.0, ancestral_state="A")
        ts1 = tables1.tree_sequence()

        tables2 = ts2.dump_tables()
        tables2.sites.add_row(position=500.0, ancestral_state="T")
        ts2 = tables2.tree_sequence()

        ta = tmc.TreesAssemblage(
            {
                tmc.ContigKey(0, 0, "chr1", "A"): ts1,
                tmc.ContigKey(1, 1, "chr2", "A"): ts2,
            }
        )
        merged = tmc.to_ts(ta)
        positions = list(merged.tables.sites.position)
        assert 500.0 in positions  # from chr1
        assert 1500.0 in positions  # from chr2 (500 + 1000)


class TestFromTreeSequence:
    def test_roundtrip(self):
        ta = make_two_contig_archive()
        ts = tmc.to_ts(ta)
        ta2 = tmc.from_ts(ts)
        assert ta2.num_contigs == 2

    def test_roundtrip_sequence_lengths(self):
        ta = make_two_contig_archive()
        ts = tmc.to_ts(ta)
        ta2 = tmc.from_ts(ts)
        assert ta2.contig("chr1").sequence_length == 1000
        assert ta2.contig("chr2").sequence_length == 2000

    def test_roundtrip_contig_metadata(self):
        ta = make_two_contig_archive()
        ts = tmc.to_ts(ta)
        ta2 = tmc.from_ts(ts)
        meta = ta2.contig("chr1").metadata[tmc.CONTIG_METADATA_KEY]
        assert meta["symbol"] == "chr1"
        assert meta["index"] == 0

    def test_no_archive_metadata_raises(self):
        tables = tskit.TableCollection(sequence_length=1000)
        tables.metadata_schema = tskit.MetadataSchema(
            {"codec": "json", "type": "object", "additionalProperties": True}
        )
        tables.metadata = {}
        ts = tables.tree_sequence()
        with pytest.raises(ValueError, match="tskit_multichrom_contigs"):
            tmc.from_ts(ts)


class TestFromSlim:
    def _make_slim_ts(self, index, id_, symbol, chrom_type, seq_len=1000):
        """Create a mock SLiM-like tree sequence."""
        tables = tskit.TableCollection(sequence_length=seq_len)
        tables.populations.add_row()
        tables.individuals.add_row()
        tables.individuals.add_row()
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0, individual=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0, individual=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0, individual=1)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0, individual=1)
        tables.nodes.add_row(flags=0, time=1.0, population=0)

        tables.edges.add_row(left=0, right=seq_len, parent=4, child=0)
        tables.edges.add_row(left=0, right=seq_len, parent=4, child=1)
        tables.edges.add_row(left=0, right=seq_len, parent=4, child=2)
        tables.edges.add_row(left=0, right=seq_len, parent=4, child=3)
        tables.sort()

        slim_meta = {
            "this_chromosome": {
                "index": index,
                "id": id_,
                "symbol": symbol,
                "type": chrom_type,
            }
        }
        schema = {
            "codec": "json",
            "type": "object",
            "additionalProperties": True,
            "properties": {
                "SLiM": {
                    "type": "object",
                    "additionalProperties": True,
                }
            },
        }
        tables.metadata_schema = tskit.MetadataSchema(schema)
        tables.metadata = {"SLiM": slim_meta}
        return tables.tree_sequence()

    def test_basic(self):
        ts1 = self._make_slim_ts(0, 0, "chr1", "A", 1000)
        ts2 = self._make_slim_ts(1, 1, "chr2", "A", 2000)
        ta = tmc.from_slim([ts1, ts2])
        assert ta.num_contigs == 2
        assert ta.contig("chr1").sequence_length == 1000
        assert ta.contig("chr2").sequence_length == 2000

    def test_nodes_marked_shared(self):
        ts1 = self._make_slim_ts(0, 0, "chr1", "A", 1000)
        ta = tmc.from_slim([ts1])
        ts = ta.contig("chr1")
        for nid in range(ts.num_nodes):
            assert ts.tables.nodes[nid].flags & tmc.NODE_IS_SHARED

    def test_contig_metadata_added(self):
        ts1 = self._make_slim_ts(0, 0, "chr1", "A", 1000)
        ta = tmc.from_slim([ts1])
        meta = ta.contig("chr1").metadata
        assert tmc.CONTIG_METADATA_KEY in meta
        assert meta[tmc.CONTIG_METADATA_KEY]["symbol"] == "chr1"

    def test_missing_slim_metadata_raises(self):
        tables = tskit.TableCollection(sequence_length=1000)
        ts = tables.tree_sequence()
        with pytest.raises(ValueError, match="SLiM"):
            tmc.from_slim([ts])

    def test_missing_this_chromosome_raises(self):
        tables = tskit.TableCollection(sequence_length=1000)
        schema = {"codec": "json", "type": "object", "additionalProperties": True}
        tables.metadata_schema = tskit.MetadataSchema(schema)
        tables.metadata = {"SLiM": {"other_key": 1}}
        ts = tables.tree_sequence()
        with pytest.raises(ValueError, match="this_chromosome"):
            tmc.from_slim([ts])
