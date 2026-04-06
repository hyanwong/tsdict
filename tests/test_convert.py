"""
Tests for conversion between TreeSequenceGroup and a single TreeSequence.
"""

import json

import pytest
import tskit

import tsgroup
from tests.conftest import make_autosomes_plus_x_archive, make_ts, make_two_contig_archive


class TestToTreeSequence:
    def test_basic(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts()
        assert isinstance(ts, tskit.TreeSequence)
        assert ts.sequence_length == tsd.total_sequence_length

    def test_metadata_present(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts()
        assert "contigs" in ts.metadata

    def test_contigs_in_metadata(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts()
        meta = ts.metadata["contigs"]
        assert len(meta) == 2
        symbols = [c["contig"]["symbol"] for c in meta]
        assert "chr1" in symbols
        assert "chr2" in symbols

    def test_empty_raises(self):
        tsd = tsgroup.TreeSequenceGroup({})
        with pytest.raises(ValueError, match="empty"):
            tsd.to_ts()

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
        tsd = tsgroup.TreeSequenceGroup(
            {
                tsgroup.ContigKey(0, 0, "c1", "A"): ts1,
                tsgroup.ContigKey(1, 1, "c2", "A"): ts2_mod,
            }
        )
        with pytest.raises(ValueError, match="[Ss]ite"):
            tsd.to_ts()

    def test_shared_nodes_have_same_id(self):
        """Shared nodes should retain their IDs in the merged tree sequence."""
        tsd = make_two_contig_archive(mark_shared=True)
        ts = tsd.to_ts()

        # Nodes 0..4 (4 samples + 1 ancestral) should all have IS_SHARED set
        for node_id in range(5):
            flags = ts.tables.nodes[node_id].flags
            assert flags & tsgroup.NODE_IS_SHARED

    def test_sequence_length(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts()
        assert ts.sequence_length == 3000  # 1000 + 2000

    def test_record_provenance_false(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts(record_provenance=False)
        assert ts.num_provenances == tsd.contig("chr1").num_provenances

    def test_record_provenance_true_adds_entry(self):
        tsd = make_two_contig_archive()
        ref_count = tsd.contig("chr1").num_provenances
        ts = tsd.to_ts(record_provenance=True)
        assert ts.num_provenances == ref_count + 1

    def test_record_provenance_true_entry_has_version(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts(record_provenance=True)
        record = json.loads(ts.tables.provenances[-1].record)
        assert record["software"]["name"] == "tsgroup"
        assert record["software"]["version"] == tsgroup.__version__

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

        tsd = tsgroup.TreeSequenceGroup(
            {
                tsgroup.ContigKey(0, 0, "chr1", "A"): ts1,
                tsgroup.ContigKey(1, 1, "chr2", "A"): ts2,
            }
        )
        merged = tsd.to_ts()
        positions = list(merged.tables.sites.position)
        assert 500.0 in positions  # from chr1
        assert 1500.0 in positions  # from chr2 (500 + 1000)

    def test_is_vacant_bits_encoded_for_mixed_sample_sets(self):
        """Verify is_vacant bits are set when contigs have different sample sets."""
        tsd = make_autosomes_plus_x_archive()

        # autosome_only_samples: shared nodes present on chr1/chr2 but not chrX
        auto_samples = set(tsd.contig("chr1").samples())
        x_samples = set(tsd.contig("chrX").samples())
        autosome_only_samples = auto_samples - x_samples
        assert len(autosome_only_samples) > 0

        # Merge and verify is_vacant bits are set for missing chrX
        merged = tsd.to_ts()
        x_index = 2  # chrX index in make_autosomes_plus_x_archive

        for nid in autosome_only_samples:
            node_meta = merged.tables.nodes[int(nid)].metadata
            assert "is_vacant" in node_meta, (
                f"Node {nid} should have is_vacant metadata; "
                f"got metadata: {node_meta}"
            )
            assert node_meta["is_vacant"] & (1 << x_index), (
                f"Node {nid} should have is_vacant bit set for chrX (index {x_index}); "
                f"got is_vacant: {node_meta['is_vacant']:032b}"
            )


class TestFromTreeSequence:
    def test_roundtrip(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts()
        tsd2 = tsgroup.from_ts(ts)
        assert tsd2.num_contigs == 2

    def test_roundtrip_sequence_lengths(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts()
        tsd2 = tsgroup.from_ts(ts)
        assert tsd2.contig("chr1").sequence_length == 1000
        assert tsd2.contig("chr2").sequence_length == 2000

    def test_roundtrip_contig_metadata(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts()
        tsd2 = tsgroup.from_ts(ts)
        meta = tsd2.contig("chr1").metadata[tsgroup.CONTIG_METADATA_KEY]
        assert meta["symbol"] == "chr1"
        assert meta["index"] == 0

    def test_record_provenance_kwarg(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts()
        tsd2 = tsgroup.from_ts(ts, record_provenance=False)
        assert tsd2.num_contigs == 2

    def test_from_ts_record_provenance_false(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts(record_provenance=False)
        tsd2 = tsgroup.from_ts(ts, record_provenance=False)
        assert tsd2.contig("chr1").num_provenances == 0
        assert tsd2.contig("chr2").num_provenances == 0

    def test_from_ts_record_provenance_true(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts(record_provenance=False)
        tsd2 = tsgroup.from_ts(ts, record_provenance=True)
        assert tsd2.contig("chr1").num_provenances == 1
        assert tsd2.contig("chr2").num_provenances == 1

    def test_from_ts_record_provenance_true_entry_has_version(self):
        tsd = make_two_contig_archive()
        ts = tsd.to_ts(record_provenance=False)
        tsd2 = tsgroup.from_ts(ts, record_provenance=True)
        record = json.loads(tsd2.contig("chr1").tables.provenances[-1].record)
        assert record["software"]["name"] == "tsgroup"
        assert record["software"]["version"] == tsgroup.__version__

    def test_no_archive_metadata_raises(self):
        tables = tskit.TableCollection(sequence_length=1000)
        tables.metadata_schema = tskit.MetadataSchema(
            {"codec": "json", "type": "object", "additionalProperties": True}
        )
        tables.metadata = {}
        ts = tables.tree_sequence()
        with pytest.raises(ValueError, match="contigs"):
            tsgroup.from_ts(ts)

    def test_roundtrip_is_vacant_bits_reconstruct_samples(self):
        """Verify is_vacant bits properly reconstruct sample markings in roundtrip."""
        # Create A/A/X archive with different sample sets per chromosome type
        tsd = make_autosomes_plus_x_archive()

        # Record original sample sets per contig
        orig_samples = {key: set(tsd[key].samples()) for key in tsd.contigs}

        # Roundtrip: to_ts encodes is_vacant bits, from_ts should decode them
        merged = tsd.to_ts()
        reconstructed = tsgroup.from_ts(merged)

        # Verify sample markings are correctly restored
        for key in tsd.contigs:
            recon_samples = set(reconstructed[key].samples())
            assert recon_samples == orig_samples[key], (
                f"Sample markings differ for {key.symbol}: "
                f"orig={orig_samples[key]}, recon={recon_samples}"
            )


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
        tsd = tsgroup.from_slim([ts1, ts2])
        assert tsd.num_contigs == 2
        assert tsd.contig("chr1").sequence_length == 1000
        assert tsd.contig("chr2").sequence_length == 2000

    def test_nodes_marked_shared(self):
        ts1 = self._make_slim_ts(0, 0, "chr1", "A", 1000)
        tsd = tsgroup.from_slim([ts1])
        ts = tsd.contig("chr1")
        for nid in range(ts.num_nodes):
            assert ts.tables.nodes[nid].flags & tsgroup.NODE_IS_SHARED

    def test_contig_metadata_added(self):
        ts1 = self._make_slim_ts(0, 0, "chr1", "A", 1000)
        tsd = tsgroup.from_slim([ts1])
        meta = tsd.contig("chr1").metadata
        assert tsgroup.CONTIG_METADATA_KEY in meta
        assert meta[tsgroup.CONTIG_METADATA_KEY]["symbol"] == "chr1"

    def test_missing_slim_metadata_raises(self):
        tables = tskit.TableCollection(sequence_length=1000)
        ts = tables.tree_sequence()
        with pytest.raises(ValueError, match="SLiM"):
            tsgroup.from_slim([ts])

    def test_missing_this_chromosome_raises(self):
        tables = tskit.TableCollection(sequence_length=1000)
        schema = {"codec": "json", "type": "object", "additionalProperties": True}
        tables.metadata_schema = tskit.MetadataSchema(schema)
        tables.metadata = {"SLiM": {"other_key": 1}}
        ts = tables.tree_sequence()
        with pytest.raises(ValueError, match="this_chromosome"):
            tsgroup.from_slim([ts])


class TestFromTreeSequences:
    def test_basic(self):
        """Test basic construction with default shared_nodes=None."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        ts2 = make_ts(seq_len=2000, num_samples=4, mark_shared=False)
        tsd = tsgroup.from_tree_sequences(
            [ts1, ts2],
            ids=[20, 21],
            symbols=["20", "21"],
            types=["A", "A"],
        )
        assert tsd.num_contigs == 2
        assert tsd.contig("20").sequence_length == 1000
        assert tsd.contig("21").sequence_length == 2000

    def test_default_shared_nodes_none(self):
        """Test that the default shared_nodes=None marks no nodes as IS_SHARED."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        tsd = tsgroup.from_tree_sequences(
            [ts1],
            ids=[1],
            symbols=["chr1"],
            types=["A"],
        )
        ts = tsd.contig("chr1")
        # Verify that no nodes have IS_SHARED marked
        for nid in range(ts.num_nodes):
            assert not (ts.tables.nodes[nid].flags & tsgroup.NODE_IS_SHARED)

    def test_contig_keys_created(self):
        """Test that ContigKey objects are created with correct values."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        ts2 = make_ts(seq_len=2000, num_samples=4, mark_shared=False)
        tsd = tsgroup.from_tree_sequences(
            [ts1, ts2],
            ids=[20, 21],
            symbols=["chr20", "chr21"],
            types=["A", "A"],
            indexes=[10, 11],
            shared_nodes="samples",
        )
        key1 = tsgroup.ContigKey(index=10, id=20, symbol="chr20", type="A")
        key2 = tsgroup.ContigKey(index=11, id=21, symbol="chr21", type="A")
        assert key1 in tsd
        assert key2 in tsd

    def test_sample_nodes_marked_shared(self):
        """Test that sample nodes are marked as IS_SHARED when shared_nodes='samples'."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        tsd = tsgroup.from_tree_sequences(
            [ts1],
            ids=[1],
            symbols=["chr1"],
            types=["A"],
            shared_nodes="samples",
        )
        ts = tsd.contig("chr1")
        sample_ids = ts.samples()
        for sample_id in sample_ids:
            assert ts.tables.nodes[sample_id].flags & tsgroup.NODE_IS_SHARED

    def test_non_sample_nodes_not_marked_shared(self):
        """Test that non-sample nodes are not marked as IS_SHARED when shared_nodes='samples'."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        tsd = tsgroup.from_tree_sequences(
            [ts1],
            ids=[1],
            symbols=["chr1"],
            types=["A"],
            shared_nodes="samples",
        )
        ts = tsd.contig("chr1")
        sample_ids = set(ts.samples())
        for nid in range(ts.num_nodes):
            if nid not in sample_ids:
                assert not (ts.tables.nodes[nid].flags & tsgroup.NODE_IS_SHARED)

    def test_list_of_node_ids_marked_shared(self):
        """Test marking a specific list of node IDs as IS_SHARED."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        # ts1 has 5 nodes: 4 samples and 1 ancestor
        shared_node_ids = [0, 1, 4]  # Mark first two samples and the ancestor
        tsd = tsgroup.from_tree_sequences(
            [ts1],
            ids=[1],
            symbols=["chr1"],
            types=["A"],
            shared_nodes=shared_node_ids,
        )
        ts = tsd.contig("chr1")
        # Check that specified nodes have IS_SHARED set
        for nid in shared_node_ids:
            assert ts.tables.nodes[nid].flags & tsgroup.NODE_IS_SHARED, f"Node {nid} not marked"
        # Check that unspecified nodes do not have IS_SHARED
        for nid in [2, 3]:
            assert not (ts.tables.nodes[nid].flags & tsgroup.NODE_IS_SHARED), f"Node {nid} marked"

    def test_default_indexes(self):
        """Test that indexes default to [0, 1, 2, ...] when not provided."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        ts2 = make_ts(seq_len=2000, num_samples=4, mark_shared=False)
        tsd = tsgroup.from_tree_sequences(
            [ts1, ts2],
            ids=[20, 21],
            symbols=["20", "21"],
            types=["A", "A"],
        )
        # Check that contigs are accessible and in order
        contigs = tsd.contigs
        assert contigs[0].index == 0
        assert contigs[1].index == 1

    def test_custom_indexes(self):
        """Test that custom indexes are respected."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        ts2 = make_ts(seq_len=2000, num_samples=4, mark_shared=False)
        tsd = tsgroup.from_tree_sequences(
            [ts1, ts2],
            ids=[20, 21],
            symbols=["20", "21"],
            types=["A", "A"],
            indexes=[100, 101],
        )
        contigs = tsd.contigs
        assert contigs[0].index == 100
        assert contigs[1].index == 101

    def test_empty_list_raises(self):
        """Test that empty tree_sequences list raises ValueError."""
        with pytest.raises(ValueError, match="empty"):
            tsgroup.from_tree_sequences(
                [],
                ids=[],
                symbols=[],
                types=[],
            )

    def test_mismatched_ids_raises(self):
        """Test that mismatched ids length raises ValueError."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        ts2 = make_ts(seq_len=2000, num_samples=4, mark_shared=False)
        with pytest.raises(ValueError, match="ids.*same length"):
            tsgroup.from_tree_sequences(
                [ts1, ts2],
                ids=[20],  # Only one id instead of two
                symbols=["20", "21"],
                types=["A", "A"],
            )

    def test_mismatched_symbols_raises(self):
        """Test that mismatched symbols length raises ValueError."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        ts2 = make_ts(seq_len=2000, num_samples=4, mark_shared=False)
        with pytest.raises(ValueError, match="symbols.*same length"):
            tsgroup.from_tree_sequences(
                [ts1, ts2],
                ids=[20, 21],
                symbols=["20"],  # Only one symbol instead of two
                types=["A", "A"],
            )

    def test_mismatched_types_raises(self):
        """Test that mismatched types length raises ValueError."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        ts2 = make_ts(seq_len=2000, num_samples=4, mark_shared=False)
        with pytest.raises(ValueError, match="types.*same length"):
            tsgroup.from_tree_sequences(
                [ts1, ts2],
                ids=[20, 21],
                symbols=["20", "21"],
                types=["A"],  # Only one type instead of two
            )

    def test_mismatched_indexes_raises(self):
        """Test that mismatched indexes length raises ValueError."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        ts2 = make_ts(seq_len=2000, num_samples=4, mark_shared=False)
        with pytest.raises(ValueError, match="indexes.*same length"):
            tsgroup.from_tree_sequences(
                [ts1, ts2],
                ids=[20, 21],
                symbols=["20", "21"],
                types=["A", "A"],
                indexes=[0],  # Only one index instead of two
            )

    def test_invalid_shared_nodes_string_raises(self):
        """Test that invalid shared_nodes string raises ValueError."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        with pytest.raises(ValueError, match="'samples'"):
            tsgroup.from_tree_sequences(
                [ts1],
                ids=[1],
                symbols=["chr1"],
                types=["A"],
                shared_nodes="invalid",
            )

    def test_invalid_shared_nodes_type_raises(self):
        """Test that non-list shared_nodes raises ValueError."""
        ts1 = make_ts(seq_len=1000, num_samples=4, mark_shared=False)
        with pytest.raises(ValueError, match="list of node IDs"):
            tsgroup.from_tree_sequences(
                [ts1],
                ids=[1],
                symbols=["chr1"],
                types=["A"],
                shared_nodes=123,  # Neither "samples" nor list-like
            )

