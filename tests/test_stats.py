"""Tests for TreeSequenceGroup stats methods."""

import numpy as np
import msprime
import pytest

import tsgroup
from tests.conftest import (
    make_autosomes_plus_x_archive,
    make_ts,
    make_two_contig_archive,
)


def _make_two_contig_archive_shared_samples_only():
    ts1 = make_ts(
        seq_len=1000,
        num_samples=4,
        contig_meta={"index": 0, "id": 0, "symbol": "chr1", "type": "A"},
        mark_shared=False,
    )
    ts2 = make_ts(
        seq_len=2000,
        num_samples=4,
        contig_meta={"index": 1, "id": 1, "symbol": "chr2", "type": "A"},
        mark_shared=False,
    )

    result = {}
    for key, ts in [
        (tsgroup.ContigKey(0, 0, "chr1", "A"), ts1),
        (tsgroup.ContigKey(1, 1, "chr2", "A"), ts2),
    ]:
        tables = ts.dump_tables()
        flags = tables.nodes.flags.copy()
        sample_ids = ts.samples()
        flags[sample_ids] = flags[sample_ids] | tsgroup.NODE_IS_SHARED
        tables.nodes.flags = flags
        result[key] = tables.tree_sequence()

    return tsgroup.TreeSequenceGroup(result)


def _add_mutations_to_archive(tsg, rate=0.05):
    """Add mutations per contig so site statistics are non-zero."""
    mutated = {}
    for i, (key, ts) in enumerate(tsg.items()):
        mutated[key] = msprime.sim_mutations(ts, rate=rate, random_seed=1000 + i)
    return tsgroup.TreeSequenceGroup(mutated)


class TestDiversity:
    @pytest.mark.parametrize("mode", ["branch", "site"])
    def test_diversity_matches_to_ts_default(self, mode):
        tsg = _make_two_contig_archive_shared_samples_only()
        tsg = _add_mutations_to_archive(tsg)

        default_sample_sets = [sorted(tsg.global_phased_node_ids)]
        expected = tsg.to_ts().diversity(sample_sets=default_sample_sets, mode=mode)
        actual = tsg.stats.diversity(mode=mode)

        assert actual == pytest.approx(expected)
        assert np.all(np.asarray(actual) > 0)

    def test_diversity_matches_to_ts_sample_sets(self):
        tsg = _make_two_contig_archive_shared_samples_only()
        sample_sets = [[0, 1], [2, 3]]
        expected = tsg.to_ts().diversity(sample_sets=sample_sets)
        actual = tsg.stats.diversity(sample_sets=sample_sets)
        assert actual == pytest.approx(expected)

    def test_diversity_matches_to_ts_span_normalise_false(self):
        tsg = _make_two_contig_archive_shared_samples_only()
        sample_sets = [sorted(tsg.global_phased_node_ids)]
        expected = tsg.to_ts().diversity(
            sample_sets=sample_sets,
            span_normalise=False,
        )
        actual = tsg.stats.diversity(
            sample_sets=sample_sets,
            span_normalise=False,
        )
        assert actual == pytest.approx(expected)

    def test_diversity_windows_not_implemented(self):
        tsg = _make_two_contig_archive_shared_samples_only()
        with pytest.raises(NotImplementedError, match="windows"):
            tsg.stats.diversity(windows=[0, 100])

    def test_diversity_empty_raises(self):
        tsg = tsgroup.TreeSequenceGroup({})
        with pytest.raises(ValueError, match="empty"):
            tsg.stats.diversity()

    def test_diversity_default_requires_no_nonglobal_samples(self):
        tsg = make_two_contig_archive(mark_shared=False)
        with pytest.raises(ValueError, match="sample_sets must be provided"):
            tsg.stats.diversity()

    def test_diversity_sample_sets_must_be_globally_phased(self):
        tsg = _make_two_contig_archive_shared_samples_only()
        with pytest.raises(ValueError, match="global_phased_node_ids"):
            tsg.stats.diversity(sample_sets=[[0, 4]])

    def test_diversity_on_autosomes_plus_x_requires_subset(self):
        tsg = make_autosomes_plus_x_archive()

        # Full A/A/X group has nonglobal sample nodes (from chrX),
        # so default sample_sets is disallowed.
        with pytest.raises(ValueError, match="sample_sets must be provided"):
            tsg.stats.diversity()

        # Restricting to autosomes removes nonglobal samples, so default works.
        tsd_auto = tsg.subset(type="A")
        value = tsd_auto.stats.diversity()
        assert value is not None

    def test_diversity_mixed_shared_samples_then_simplify_matches(self):
        # Mixed-sample group: some sample IDs are globally shared,
        # others are not (partial-sample ARG across contigs).
        tsg = make_autosomes_plus_x_archive()
        tsg = _add_mutations_to_archive(tsg)

        with pytest.raises(ValueError, match="sample_sets must be provided"):
            tsg.stats.diversity()

        shared_samples = sorted(
            s for s in tsg.global_phased_node_ids if tsg.contig("chr1").node(s).is_sample()
        )
        explicit_value = tsg.stats.diversity(sample_sets=[shared_samples])

        tsd_shared = tsg.simplify(samples=shared_samples)
        assert not tsd_shared.is_nonglobal_sample_arg
        simplified_value = tsd_shared.stats.diversity()

        assert simplified_value == pytest.approx(explicit_value)


class TestPCA:
    def test_pca_default_matches_to_ts(self):
        tsg = _make_two_contig_archive_shared_samples_only()
        samples = np.asarray(sorted(tsg.global_phased_node_ids), dtype=np.int32)
        expected = tsg.to_ts().pca(num_components=2, samples=samples, random_seed=42)
        actual = tsg.stats.pca(num_components=2, random_seed=42)

        assert actual.eigenvalues == pytest.approx(expected.eigenvalues)
        assert np.abs(actual.factors) == pytest.approx(np.abs(expected.factors))

    def test_pca_with_explicit_samples(self):
        tsg = _make_two_contig_archive_shared_samples_only()
        samples = sorted(tsg.global_phased_node_ids)[:2]
        result = tsg.stats.pca(num_components=1, samples=samples, random_seed=1)
        assert result.factors.shape == (2, 1)

    def test_pca_samples_must_be_globally_phased(self):
        tsg = make_autosomes_plus_x_archive()
        nonglobal = sorted(
            set(range(tsg.contig("chr1").num_samples)) - tsg.global_phased_node_ids
        )
        assert len(nonglobal) > 0, "expected nonglobal sample nodes"
        with pytest.raises(ValueError, match="global_phased_node_ids"):
            tsg.stats.pca(num_components=1, samples=nonglobal, random_seed=1)

    def test_pca_default_raises_for_nonglobal_sample_arg(self):
        tsg = make_autosomes_plus_x_archive()
        assert tsg.is_nonglobal_sample_arg
        with pytest.raises(ValueError, match="individuals or samples"):
            tsg.stats.pca(num_components=1)

    def test_pca_with_individuals_rejects_mixed_types(self):
        tsg = make_autosomes_plus_x_archive()
        assert tsg.is_nonglobal_sample_arg
        with pytest.raises(ValueError, match="same type"):
            tsg.stats.pca(num_components=1, individuals=[0, 1])

    def test_pca_with_individuals_matches_to_ts(self):
        tsg = _make_two_contig_archive_shared_samples_only()
        individuals = [0, 1]
        expected = tsg.to_ts().pca(
            num_components=1, individuals=individuals, random_seed=42
        )
        actual = tsg.stats.pca(
            num_components=1, individuals=individuals, random_seed=42
        )
        assert actual.eigenvalues == pytest.approx(expected.eigenvalues)
        assert np.abs(actual.factors) == pytest.approx(np.abs(expected.factors))

    def test_pca_cannot_specify_samples_and_individuals(self):
        tsg = _make_two_contig_archive_shared_samples_only()
        with pytest.raises(ValueError, match="Cannot specify both"):
            tsg.stats.pca(num_components=1, samples=[0], individuals=[0])

    def test_pca_with_individuals_fully_nonglobal_sample_arg(self):
        tsg = make_two_contig_archive(mark_shared=False)
        assert tsg.is_nonglobal_sample_arg
        assert len(tsg.global_phased_node_ids) == 0

        with pytest.raises(ValueError, match="individuals or samples"):
            tsg.stats.pca(num_components=1)

        individuals = [0, 1]
        ts_compact = tsg.to_ts()
        ts_compact = ts_compact.simplify(
            samples=ts_compact.samples(), record_provenance=False
        )
        expected = ts_compact.pca(
            num_components=1, individuals=individuals, random_seed=7
        )
        actual = tsg.stats.pca(
            num_components=1, individuals=individuals, random_seed=7
        )
        assert actual.eigenvalues == pytest.approx(expected.eigenvalues)
        assert np.abs(actual.factors) == pytest.approx(np.abs(expected.factors))
