"""Tests for TreesAssemblage stats methods."""

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

    return tmc.TreesAssemblage(result)


def _add_mutations_to_archive(ta, rate=0.05):
    """Add mutations per contig so site statistics are non-zero."""
    mutated = {}
    for i, (key, ts) in enumerate(ta.items()):
        mutated[key] = msprime.sim_mutations(ts, rate=rate, random_seed=1000 + i)
    return tmc.TreesAssemblage(mutated)


class TestDiversity:
    @pytest.mark.parametrize("mode", ["branch", "site"])
    def test_diversity_matches_to_ts_default(self, mode):
        ta = _make_two_contig_archive_shared_samples_only()
        ta = _add_mutations_to_archive(ta)

        default_sample_sets = [sorted(ta.global_phased_node_ids)]
        expected = ta.to_ts().diversity(sample_sets=default_sample_sets, mode=mode)
        actual = ta.stats.diversity(mode=mode)

        assert actual == pytest.approx(expected)
        assert np.all(np.asarray(actual) > 0)

    def test_diversity_matches_to_ts_sample_sets(self):
        ta = _make_two_contig_archive_shared_samples_only()
        sample_sets = [[0, 1], [2, 3]]
        expected = ta.to_ts().diversity(sample_sets=sample_sets)
        actual = ta.stats.diversity(sample_sets=sample_sets)
        assert actual == pytest.approx(expected)

    def test_diversity_matches_to_ts_span_normalise_false(self):
        ta = _make_two_contig_archive_shared_samples_only()
        sample_sets = [sorted(ta.global_phased_node_ids)]
        expected = ta.to_ts().diversity(
            sample_sets=sample_sets,
            span_normalise=False,
        )
        actual = ta.stats.diversity(
            sample_sets=sample_sets,
            span_normalise=False,
        )
        assert actual == pytest.approx(expected)

    def test_diversity_windows_not_implemented(self):
        ta = _make_two_contig_archive_shared_samples_only()
        with pytest.raises(NotImplementedError, match="windows"):
            ta.stats.diversity(windows=[0, 100])

    def test_diversity_empty_raises(self):
        ta = tmc.TreesAssemblage({})
        with pytest.raises(ValueError, match="empty"):
            ta.stats.diversity()

    def test_diversity_default_requires_no_nonglobal_samples(self):
        ta = make_two_contig_archive(mark_shared=False)
        with pytest.raises(ValueError, match="sample_sets must be provided"):
            ta.stats.diversity()

    def test_diversity_sample_sets_must_be_globally_phased(self):
        ta = _make_two_contig_archive_shared_samples_only()
        with pytest.raises(ValueError, match="global_phased_node_ids"):
            ta.stats.diversity(sample_sets=[[0, 4]])

    def test_diversity_on_autosomes_plus_x_requires_subset(self):
        ta = make_autosomes_plus_x_archive()

        # Full A/A/X assemblage has nonglobal sample nodes (from chrX),
        # so default sample_sets is disallowed.
        with pytest.raises(ValueError, match="sample_sets must be provided"):
            ta.stats.diversity()

        # Restricting to autosomes removes nonglobal samples, so default works.
        ta_auto = ta.subset(type="A")
        value = ta_auto.stats.diversity()
        assert value is not None

    def test_diversity_mixed_shared_samples_then_simplify_matches(self):
        # Mixed-sample assemblage: some sample IDs are globally shared,
        # others are not (partial-sample ARG across contigs).
        ta = make_autosomes_plus_x_archive()
        ta = _add_mutations_to_archive(ta)

        with pytest.raises(ValueError, match="sample_sets must be provided"):
            ta.stats.diversity()

        shared_samples = sorted(
            s for s in ta.global_phased_node_ids if ta.contig("chr1").node(s).is_sample()
        )
        explicit_value = ta.stats.diversity(sample_sets=[shared_samples])

        ta_shared = ta.simplify(samples=shared_samples)
        assert not ta_shared.is_nonglobal_sample_arg
        simplified_value = ta_shared.stats.diversity()

        assert simplified_value == pytest.approx(explicit_value)
