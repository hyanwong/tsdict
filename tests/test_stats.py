"""Tests for TreeSequenceDictionary stats methods."""

import numpy as np
import msprime
import pytest

import tskit_multichrom as tmc
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
        (tmc.ContigKey(0, 0, "chr1", "A"), ts1),
        (tmc.ContigKey(1, 1, "chr2", "A"), ts2),
    ]:
        tables = ts.dump_tables()
        flags = tables.nodes.flags.copy()
        sample_ids = ts.samples()
        flags[sample_ids] = flags[sample_ids] | tmc.NODE_IS_SHARED
        tables.nodes.flags = flags
        result[key] = tables.tree_sequence()

    return tmc.TreeSequenceDictionary(result)


def _add_mutations_to_archive(tsd, rate=0.05):
    """Add mutations per contig so site statistics are non-zero."""
    mutated = {}
    for i, (key, ts) in enumerate(tsd.items()):
        mutated[key] = msprime.sim_mutations(ts, rate=rate, random_seed=1000 + i)
    return tmc.TreeSequenceDictionary(mutated)


class TestDiversity:
    @pytest.mark.parametrize("mode", ["branch", "site"])
    def test_diversity_matches_to_ts_default(self, mode):
        tsd = _make_two_contig_archive_shared_samples_only()
        tsd = _add_mutations_to_archive(tsd)

        default_sample_sets = [sorted(tsd.global_phased_node_ids)]
        expected = tsd.to_ts().diversity(sample_sets=default_sample_sets, mode=mode)
        actual = tsd.stats.diversity(mode=mode)

        assert actual == pytest.approx(expected)
        assert np.all(np.asarray(actual) > 0)

    def test_diversity_matches_to_ts_sample_sets(self):
        tsd = _make_two_contig_archive_shared_samples_only()
        sample_sets = [[0, 1], [2, 3]]
        expected = tsd.to_ts().diversity(sample_sets=sample_sets)
        actual = tsd.stats.diversity(sample_sets=sample_sets)
        assert actual == pytest.approx(expected)

    def test_diversity_matches_to_ts_span_normalise_false(self):
        tsd = _make_two_contig_archive_shared_samples_only()
        sample_sets = [sorted(tsd.global_phased_node_ids)]
        expected = tsd.to_ts().diversity(
            sample_sets=sample_sets,
            span_normalise=False,
        )
        actual = tsd.stats.diversity(
            sample_sets=sample_sets,
            span_normalise=False,
        )
        assert actual == pytest.approx(expected)

    def test_diversity_windows_not_implemented(self):
        tsd = _make_two_contig_archive_shared_samples_only()
        with pytest.raises(NotImplementedError, match="windows"):
            tsd.stats.diversity(windows=[0, 100])

    def test_diversity_empty_raises(self):
        tsd = tmc.TreeSequenceDictionary({})
        with pytest.raises(ValueError, match="empty"):
            tsd.stats.diversity()

    def test_diversity_default_requires_no_nonglobal_samples(self):
        tsd = make_two_contig_archive(mark_shared=False)
        with pytest.raises(ValueError, match="sample_sets must be provided"):
            tsd.stats.diversity()

    def test_diversity_sample_sets_must_be_globally_phased(self):
        tsd = _make_two_contig_archive_shared_samples_only()
        with pytest.raises(ValueError, match="global_phased_node_ids"):
            tsd.stats.diversity(sample_sets=[[0, 4]])

    def test_diversity_on_autosomes_plus_x_requires_subset(self):
        tsd = make_autosomes_plus_x_archive()

        # Full A/A/X assemblage has nonglobal sample nodes (from chrX),
        # so default sample_sets is disallowed.
        with pytest.raises(ValueError, match="sample_sets must be provided"):
            tsd.stats.diversity()

        # Restricting to autosomes removes nonglobal samples, so default works.
        tsd_auto = tsd.subset(type="A")
        value = tsd_auto.stats.diversity()
        assert value is not None

    def test_diversity_mixed_shared_samples_then_simplify_matches(self):
        # Mixed-sample assemblage: some sample IDs are globally shared,
        # others are not (partial-sample ARG across contigs).
        tsd = make_autosomes_plus_x_archive()
        tsd = _add_mutations_to_archive(tsd)

        with pytest.raises(ValueError, match="sample_sets must be provided"):
            tsd.stats.diversity()

        shared_samples = sorted(
            s for s in tsd.global_phased_node_ids if tsd.contig("chr1").node(s).is_sample()
        )
        explicit_value = tsd.stats.diversity(sample_sets=[shared_samples])

        tsd_shared = tsd.simplify(samples=shared_samples)
        assert not tsd_shared.is_nonglobal_sample_arg
        simplified_value = tsd_shared.stats.diversity()

        assert simplified_value == pytest.approx(explicit_value)


class TestPCA:
    def test_pca_default_matches_to_ts(self):
        tsd = _make_two_contig_archive_shared_samples_only()
        samples = np.asarray(sorted(tsd.global_phased_node_ids), dtype=np.int32)
        expected = tsd.to_ts().pca(num_components=2, samples=samples, random_seed=42)
        actual = tsd.stats.pca(num_components=2, random_seed=42)

        assert actual.eigenvalues == pytest.approx(expected.eigenvalues)
        assert np.abs(actual.factors) == pytest.approx(np.abs(expected.factors))

    def test_pca_with_explicit_samples(self):
        tsd = _make_two_contig_archive_shared_samples_only()
        samples = sorted(tsd.global_phased_node_ids)[:2]
        result = tsd.stats.pca(num_components=1, samples=samples, random_seed=1)
        assert result.factors.shape == (2, 1)

    def test_pca_samples_must_be_globally_phased(self):
        tsd = make_autosomes_plus_x_archive()
        nonglobal = sorted(
            set(range(tsd.contig("chr1").num_samples)) - tsd.global_phased_node_ids
        )
        assert len(nonglobal) > 0, "expected nonglobal sample nodes"
        with pytest.raises(ValueError, match="global_phased_node_ids"):
            tsd.stats.pca(num_components=1, samples=nonglobal, random_seed=1)

    def test_pca_default_raises_for_nonglobal_sample_arg(self):
        tsd = make_autosomes_plus_x_archive()
        assert tsd.is_nonglobal_sample_arg
        with pytest.raises(ValueError, match="individuals or samples"):
            tsd.stats.pca(num_components=1)

    def test_pca_with_individuals_rejects_mixed_types(self):
        tsd = make_autosomes_plus_x_archive()
        assert tsd.is_nonglobal_sample_arg
        with pytest.raises(ValueError, match="same type"):
            tsd.stats.pca(num_components=1, individuals=[0, 1])

    def test_pca_with_individuals_matches_to_ts(self):
        tsd = _make_two_contig_archive_shared_samples_only()
        individuals = [0, 1]
        expected = tsd.to_ts().pca(
            num_components=1, individuals=individuals, random_seed=42
        )
        actual = tsd.stats.pca(
            num_components=1, individuals=individuals, random_seed=42
        )
        assert actual.eigenvalues == pytest.approx(expected.eigenvalues)
        assert np.abs(actual.factors) == pytest.approx(np.abs(expected.factors))

    def test_pca_cannot_specify_samples_and_individuals(self):
        tsd = _make_two_contig_archive_shared_samples_only()
        with pytest.raises(ValueError, match="Cannot specify both"):
            tsd.stats.pca(num_components=1, samples=[0], individuals=[0])

    def test_pca_with_individuals_fully_nonglobal_sample_arg(self):
        tsd = make_two_contig_archive(mark_shared=False)
        assert tsd.is_nonglobal_sample_arg
        assert len(tsd.global_phased_node_ids) == 0

        with pytest.raises(ValueError, match="individuals or samples"):
            tsd.stats.pca(num_components=1)

        individuals = [0, 1]
        ts_compact = tsd.to_ts()
        ts_compact = ts_compact.simplify(
            samples=ts_compact.samples(), record_provenance=False
        )
        expected = ts_compact.pca(
            num_components=1, individuals=individuals, random_seed=7
        )
        actual = tsd.stats.pca(
            num_components=1, individuals=individuals, random_seed=7
        )
        assert actual.eigenvalues == pytest.approx(expected.eigenvalues)
        assert np.abs(actual.factors) == pytest.approx(np.abs(expected.factors))
