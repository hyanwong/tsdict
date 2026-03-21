"""
Tests for the TreesArchive core class.
"""

import pytest
import tskit

import tskit_multichrom as tmc
from tests.conftest import make_ts, make_two_contig_archive


class TestContigKey:
    def test_basic(self):
        key = tmc.ContigKey(0, 0, "chr1", "A")
        assert key.index == 0
        assert key.id == 0
        assert key.symbol == "chr1"
        assert key.type == "A"

    def test_namedtuple(self):
        key = tmc.ContigKey(0, 0, "chr1", "A")
        assert key[0] == 0
        assert key[1] == 0
        assert key[2] == "chr1"
        assert key[3] == "A"


class TestTreesArchiveConstruction:
    def test_basic(self):
        ta = make_two_contig_archive()
        assert ta.num_contigs == 2

    def test_empty(self):
        ta = tmc.TreesArchive({})
        assert ta.num_contigs == 0

    def test_wrong_key_type(self):
        ts = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        with pytest.raises(TypeError, match="ContigKey"):
            tmc.TreesArchive({"bad_key": ts})

    def test_wrong_value_type(self):
        with pytest.raises(TypeError, match="tskit.TreeSequence"):
            tmc.TreesArchive({tmc.ContigKey(0, 0, "c1", "A"): "not a ts"})

    def test_non_dict(self):
        with pytest.raises(TypeError, match="dict"):
            tmc.TreesArchive([])


class TestValidation:
    def test_duplicate_index(self):
        ts1 = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        ts2 = make_ts(contig_meta={"index": 0, "id": 1, "symbol": "c2", "type": "A"})
        with pytest.raises(ValueError, match="index"):
            tmc.TreesArchive(
                {
                    tmc.ContigKey(0, 0, "c1", "A"): ts1,
                    tmc.ContigKey(0, 1, "c2", "A"): ts2,
                }
            )

    def test_duplicate_id(self):
        ts1 = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        ts2 = make_ts(contig_meta={"index": 1, "id": 0, "symbol": "c2", "type": "A"})
        with pytest.raises(ValueError, match="id"):
            tmc.TreesArchive(
                {
                    tmc.ContigKey(0, 0, "c1", "A"): ts1,
                    tmc.ContigKey(1, 0, "c2", "A"): ts2,
                }
            )

    def test_duplicate_symbol(self):
        ts1 = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        ts2 = make_ts(contig_meta={"index": 1, "id": 1, "symbol": "c1", "type": "A"})
        with pytest.raises(ValueError, match="symbol"):
            tmc.TreesArchive(
                {
                    tmc.ContigKey(0, 0, "c1", "A"): ts1,
                    tmc.ContigKey(1, 1, "c1", "A"): ts2,
                }
            )

    def test_missing_contig_metadata(self):
        ts_no_meta = make_ts(contig_meta=None)  # no metadata at all
        with pytest.raises(ValueError, match="contig"):
            tmc.TreesArchive({tmc.ContigKey(0, 0, "c1", "A"): ts_no_meta})

    def test_individual_tables_mismatch(self):
        ts1 = make_ts(
            num_samples=4,
            contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"},
        )
        ts2 = make_ts(
            num_samples=6,  # different number of individuals
            contig_meta={"index": 1, "id": 1, "symbol": "c2", "type": "A"},
        )
        with pytest.raises(ValueError, match="Individual"):
            tmc.TreesArchive(
                {
                    tmc.ContigKey(0, 0, "c1", "A"): ts1,
                    tmc.ContigKey(1, 1, "c2", "A"): ts2,
                }
            )

    def test_population_tables_mismatch(self):
        ts1 = make_ts(
            contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"},
        )
        # Create ts with different populations
        tables = ts1.dump_tables()
        tables.populations.add_row()
        tables.metadata = {"contig": {"index": 1, "id": 1, "symbol": "c2", "type": "A"}}
        ts2 = tables.tree_sequence()
        with pytest.raises(ValueError, match="Population"):
            tmc.TreesArchive(
                {
                    tmc.ContigKey(0, 0, "c1", "A"): ts1,
                    tmc.ContigKey(1, 1, "c2", "A"): ts2,
                }
            )

    def test_with_migrations_fails(self):
        ts = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        tables = ts.dump_tables()
        # Add two populations and a migration between them
        tables.populations.add_row()
        tables.migrations.add_row(
            left=0, right=1000, node=0, source=0, dest=1, time=0.5
        )
        ts_with_mig = tables.tree_sequence()
        with pytest.raises(ValueError, match="[Mm]igration"):
            tmc.TreesArchive({tmc.ContigKey(0, 0, "c1", "A"): ts_with_mig})


class TestCaching:
    def test_total_sequence_length(self):
        ta = make_two_contig_archive()
        assert ta.total_sequence_length == 1000 + 2000

    def test_cross_phased_nodes(self):
        ta = make_two_contig_archive(mark_shared=True)
        # All nodes have IS_SHARED set, so all should be cross-phased
        assert len(ta.cross_phased_node_ids) == 5  # 4 samples + 1 ancestral

    def test_no_cross_phased_nodes(self):
        ta = make_two_contig_archive(mark_shared=False)
        assert len(ta.cross_phased_node_ids) == 0

    def test_not_partial_sample_arg_when_shared(self):
        ta = make_two_contig_archive(mark_shared=True)
        assert not ta.is_partial_sample_arg

    def test_is_partial_sample_arg_when_not_shared(self):
        ta = make_two_contig_archive(mark_shared=False)
        assert ta.is_partial_sample_arg

    def test_nonglobal_sample_count(self):
        ta = make_two_contig_archive(mark_shared=False)
        # 4 samples in each contig but they don't overlap
        assert ta.nonglobal_sample_node_count == 4


class TestAccess:
    def test_chr_access(self):
        ta = make_two_contig_archive()
        ts = ta.chr("chr1")
        assert ts.sequence_length == 1000

    def test_chr_unknown(self):
        ta = make_two_contig_archive()
        with pytest.raises(KeyError, match="chr99"):
            ta.chr("chr99")

    def test_contig_by_id(self):
        ta = make_two_contig_archive()
        ts = ta.contig(1)
        assert ts.sequence_length == 2000

    def test_contig_by_index(self):
        ta = make_two_contig_archive()
        ts = ta.contig(0)
        assert ts.sequence_length == 1000

    def test_getitem(self):
        ta = make_two_contig_archive()
        key = tmc.ContigKey(0, 0, "chr1", "A")
        assert ta[key].sequence_length == 1000

    def test_len(self):
        ta = make_two_contig_archive()
        assert len(ta) == 2

    def test_iter(self):
        ta = make_two_contig_archive()
        keys = list(ta)
        assert keys[0].symbol == "chr1"
        assert keys[1].symbol == "chr2"

    def test_repr(self):
        ta = make_two_contig_archive()
        r = repr(ta)
        assert "chr1" in r
        assert "chr2" in r


class TestSubset:
    def test_subset_by_symbol(self):
        ta = make_two_contig_archive()
        sub = ta.subset(symbols=["chr1"])
        assert sub.num_contigs == 1
        assert sub.chr("chr1").sequence_length == 1000

    def test_subset_by_type(self):
        # Create an archive with different types
        ts1 = make_ts(
            seq_len=1000,
            contig_meta={"index": 0, "id": 0, "symbol": "chr1", "type": "A"},
        )
        ts2 = make_ts(
            seq_len=500,
            contig_meta={"index": 1, "id": 1, "symbol": "chrX", "type": "X"},
        )
        ta = tmc.TreesArchive(
            {
                tmc.ContigKey(0, 0, "chr1", "A"): ts1,
                tmc.ContigKey(1, 1, "chrX", "X"): ts2,
            }
        )
        sub = ta.subset(types=["A"])
        assert sub.num_contigs == 1
        assert sub.chr("chr1")

    def test_subset_empty_result(self):
        ta = make_two_contig_archive()
        sub = ta.subset(symbols=["nonexistent"])
        assert sub.num_contigs == 0

    def test_subset_by_id(self):
        ta = make_two_contig_archive()
        sub = ta.subset(ids=[0])
        assert sub.num_contigs == 1


class TestReindex:
    def test_reindex_default(self):
        ta = make_two_contig_archive()
        ta2 = ta.reindex()
        keys = ta2.contigs
        assert [k.index for k in keys] == [0, 1]
        assert [k.symbol for k in keys] == ["chr1", "chr2"]

    def test_reindex_custom_order(self):
        ta = make_two_contig_archive()
        ta2 = ta.reindex(order=["chr2", "chr1"])
        keys = ta2.contigs
        assert keys[0].symbol == "chr2"
        assert keys[0].index == 0
        assert keys[1].symbol == "chr1"
        assert keys[1].index == 1

    def test_reindex_updates_metadata(self):
        ta = make_two_contig_archive()
        ta2 = ta.reindex(order=["chr2", "chr1"])
        meta0 = ta2.chr("chr2").metadata["contig"]
        assert meta0["index"] == 0
        meta1 = ta2.chr("chr1").metadata["contig"]
        assert meta1["index"] == 1
