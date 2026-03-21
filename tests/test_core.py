"""
Tests for the TreesAssemblage core class.
"""

import pytest
import tskit

import tskit_multichrom as tmc
from tests.conftest import make_ts, make_two_contig_archive, make_autosomes_plus_x_archive


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


class TestTreesAssemblageConstruction:
    def test_basic(self):
        ta = make_two_contig_archive()
        assert ta.num_contigs == 2

    def test_empty(self):
        ta = tmc.TreesAssemblage({})
        assert ta.num_contigs == 0

    def test_wrong_key_type(self):
        ts = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        with pytest.raises(TypeError, match="ContigKey"):
            tmc.TreesAssemblage({"bad_key": ts})

    def test_wrong_value_type(self):
        with pytest.raises(TypeError, match="tskit.TreeSequence"):
            tmc.TreesAssemblage({tmc.ContigKey(0, 0, "c1", "A"): "not a ts"})

    def test_non_dict(self):
        with pytest.raises(TypeError, match="dict"):
            tmc.TreesAssemblage([])


class TestValidation:
    def test_static_validate_accepts_mapping(self):
        ts1 = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        ts2 = make_ts(contig_meta={"index": 1, "id": 1, "symbol": "c2", "type": "A"})
        tmc.TreesAssemblage.validate(
            {
                tmc.ContigKey(0, 0, "c1", "A"): ts1,
                tmc.ContigKey(1, 1, "c2", "A"): ts2,
            }
        )

    def test_duplicate_index(self):
        ts1 = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        ts2 = make_ts(contig_meta={"index": 0, "id": 1, "symbol": "c2", "type": "A"})
        with pytest.raises(ValueError, match="index"):
            tmc.TreesAssemblage(
                {
                    tmc.ContigKey(0, 0, "c1", "A"): ts1,
                    tmc.ContigKey(0, 1, "c2", "A"): ts2,
                }
            )

    def test_duplicate_id(self):
        ts1 = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        ts2 = make_ts(contig_meta={"index": 1, "id": 0, "symbol": "c2", "type": "A"})
        with pytest.raises(ValueError, match="id"):
            tmc.TreesAssemblage(
                {
                    tmc.ContigKey(0, 0, "c1", "A"): ts1,
                    tmc.ContigKey(1, 0, "c2", "A"): ts2,
                }
            )

    def test_duplicate_symbol(self):
        ts1 = make_ts(contig_meta={"index": 0, "id": 0, "symbol": "c1", "type": "A"})
        ts2 = make_ts(contig_meta={"index": 1, "id": 1, "symbol": "c1", "type": "A"})
        with pytest.raises(ValueError, match="symbol"):
            tmc.TreesAssemblage(
                {
                    tmc.ContigKey(0, 0, "c1", "A"): ts1,
                    tmc.ContigKey(1, 1, "c1", "A"): ts2,
                }
            )

    def test_missing_contig_metadata(self):
        ts_no_meta = make_ts(contig_meta=None)  # no metadata at all
        with pytest.raises(ValueError, match="contig"):
            tmc.TreesAssemblage({tmc.ContigKey(0, 0, "c1", "A"): ts_no_meta})

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
            tmc.TreesAssemblage(
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
            tmc.TreesAssemblage(
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
            tmc.TreesAssemblage({tmc.ContigKey(0, 0, "c1", "A"): ts_with_mig})


class TestCaching:
    def test_total_sequence_length(self):
        ta = make_two_contig_archive()
        assert ta.total_sequence_length == 1000 + 2000

    def test_shared_nodes(self):
        ta = make_two_contig_archive(mark_shared=True)
        assert len(ta.shared_node_ids) == 5  # 4 samples + 1 ancestral

    def test_no_shared_nodes(self):
        ta = make_two_contig_archive(mark_shared=False)
        assert len(ta.shared_node_ids) == 0

    def test_global_phased_nodes(self):
        ta = make_two_contig_archive(mark_shared=True)
        # All nodes have IS_SHARED set, so all should be globally phased
        assert len(ta.global_phased_node_ids) == 5  # 4 samples + 1 ancestral

    def test_no_global_phased_nodes(self):
        ta = make_two_contig_archive(mark_shared=False)
        assert len(ta.global_phased_node_ids) == 0

    def test_not_nonglobal_sample_arg_when_shared(self):
        ta = make_two_contig_archive(mark_shared=True)
        assert not ta.is_nonglobal_sample_arg

    def test_is_nonglobal_sample_arg_when_not_shared(self):
        ta = make_two_contig_archive(mark_shared=False)
        assert ta.is_nonglobal_sample_arg

    def test_nonglobal_sample_count(self):
        ta = make_two_contig_archive(mark_shared=False)
        # 4 samples in each contig but they don't overlap
        assert ta.nonglobal_sample_node_count == 4


class TestAccess:
    def test_contig_by_symbol(self):
        ta = make_two_contig_archive()
        ts = ta.contig("chr1")
        assert ts.sequence_length == 1000

    def test_contig_symbol_unknown(self):
        ta = make_two_contig_archive()
        with pytest.raises(KeyError, match="chr99"):
            ta.contig("chr99")

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

    def test_keys(self):
        ta = make_two_contig_archive()
        k = ta.keys()
        assert [x.symbol for x in k] == ["chr1", "chr2"]

    def test_values(self):
        ta = make_two_contig_archive()
        v = ta.values()
        assert [ts.sequence_length for ts in v] == [1000, 2000]

    def test_items(self):
        ta = make_two_contig_archive()
        items = ta.items()
        assert len(items) == 2
        keys, tss = zip(*items)
        assert [k.symbol for k in keys] == ["chr1", "chr2"]
        assert [ts.sequence_length for ts in tss] == [1000, 2000]


class TestSubset:
    def test_subset_by_symbol(self):
        ta = make_two_contig_archive()
        sub = ta.subset(symbols=["chr1"])
        assert sub.num_contigs == 1
        assert sub.contig("chr1").sequence_length == 1000

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
        ta = tmc.TreesAssemblage(
            {
                tmc.ContigKey(0, 0, "chr1", "A"): ts1,
                tmc.ContigKey(1, 1, "chrX", "X"): ts2,
            }
        )
        sub = ta.subset(types=["A"])
        assert sub.num_contigs == 1
        assert sub.contig("chr1")

    def test_subset_by_single_type(self):
        """subset(type=...) should work as a convenient single-value filter."""
        ts1 = make_ts(
            seq_len=1000,
            contig_meta={"index": 0, "id": 0, "symbol": "chr1", "type": "A"},
        )
        ts2 = make_ts(
            seq_len=500,
            contig_meta={"index": 1, "id": 1, "symbol": "chrX", "type": "X"},
        )
        ta = tmc.TreesAssemblage(
            {
                tmc.ContigKey(0, 0, "chr1", "A"): ts1,
                tmc.ContigKey(1, 1, "chrX", "X"): ts2,
            }
        )
        sub = ta.subset(type="A")
        assert sub.num_contigs == 1
        assert sub.contig("chr1")

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
        meta0 = ta2.contig("chr2").metadata["contig"]
        assert meta0["index"] == 0
        meta1 = ta2.contig("chr1").metadata["contig"]
        assert meta1["index"] == 1


class TestSimplify:
    def test_simplify_reduces_nodes(self):
        ta = make_two_contig_archive(mark_shared=True)
        ta2 = ta.simplify()
        # Ancestral nodes should be simplified away (star topology → one root)
        for ts in ta2.values():
            assert ts.num_nodes <= ta.contig("chr1").num_nodes

    def test_simplify_preserves_num_samples(self):
        ta = make_two_contig_archive(mark_shared=True, num_samples=4)
        ta2 = ta.simplify()
        for ts in ta2.values():
            assert ts.num_samples == 4

    def test_simplify_preserves_is_shared_flags(self):
        ta = make_two_contig_archive(mark_shared=True)
        ta2 = ta.simplify()
        assert not ta2.is_nonglobal_sample_arg
        # All sample nodes should still have IS_SHARED set after simplify
        ts = ta2.contig("chr1")
        for s in ts.samples():
            assert ts.node(s).flags & tmc.NODE_IS_SHARED

    def test_simplify_consistent_sample_ids_across_contigs(self):
        ta = make_two_contig_archive(mark_shared=True, num_samples=4)
        ta2 = ta.simplify()
        # Sample node IDs must be the same in both contigs
        samples_chr1 = set(ta2.contig("chr1").samples())
        samples_chr2 = set(ta2.contig("chr2").samples())
        assert samples_chr1 == samples_chr2

    def test_simplify_fails_with_nonglobal_samples(self):
        ta = make_two_contig_archive(mark_shared=False)
        with pytest.raises(ValueError, match="nonglobal"):
            ta.simplify()

    def test_simplify_empty_assemblage(self):
        ta = tmc.TreesAssemblage({})
        ta2 = ta.simplify()
        assert ta2.num_contigs == 0

    def test_simplify_by_samples(self):
        ta = make_two_contig_archive(mark_shared=True, num_samples=4)
        global_samples = sorted(ta.global_phased_node_ids & set(ta.contig("chr1").samples()))
        ta2 = ta.simplify(samples=global_samples[:2])
        for ts in ta2.values():
            assert ts.num_samples == 2

    def test_simplify_by_individuals(self):
        # Simplify to a subset of individuals, works even without cross-phasing
        ta = make_two_contig_archive(mark_shared=True, num_samples=4)
        ta2 = ta.simplify(individuals=[0, 1])  # keep 2 out of 2 individuals
        for ts in ta2.values():
            assert ts.num_samples == 4  # 2 individuals × 2 haplotypes

    def test_simplify_by_individuals_subset(self):
        ta = make_two_contig_archive(mark_shared=True, num_samples=4)
        ta2 = ta.simplify(individuals=[0])  # keep 1 out of 2 individuals
        for ts in ta2.values():
            assert ts.num_samples == 2  # 1 individual × 2 haplotypes

    def test_simplify_by_individuals_on_partial_arg(self):
        # Can simplify a partial-sample ARG using individuals=
        ta = make_autosomes_plus_x_archive(num_samples=4)
        assert ta.is_nonglobal_sample_arg
        ta2 = ta.simplify(individuals=[0])
        assert ta2.num_contigs == 3

    def test_simplify_samples_and_individuals_raises(self):
        ta = make_two_contig_archive(mark_shared=True)
        with pytest.raises(ValueError, match="both"):
            ta.simplify(samples=[0, 1], individuals=[0])
