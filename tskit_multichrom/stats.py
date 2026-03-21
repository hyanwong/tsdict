"""Statistics helpers for TreesAssemblage."""

import numpy as np


class TreesAssemblageStats:
    """Statistics namespace for :class:`tskit_multichrom.core.TreesAssemblage`."""

    def __init__(self, ta):
        self._ta = ta

    def _resolve_sample_sets(self, sample_sets):
        """
        Resolve and validate ``sample_sets`` for multi-contig statistics.

        Rules:
        - If ``sample_sets is None`` and there are no nonglobal samples,
          default to one sample set containing all ``global_phased_node_ids``.
        - If ``sample_sets is None`` and nonglobal samples exist, raise.
        - If ``sample_sets`` is provided, every sample ID must be globally phased.
        """
        global_phased = set(self._ta.global_phased_node_ids)

        if sample_sets is None:
            if self._ta.nonglobal_sample_node_count != 0:
                raise ValueError(
                    "sample_sets must be provided when nonglobal sample nodes are present"
                )
            return [sorted(global_phased)]

        for sid in _iter_sample_ids(sample_sets):
            if sid not in global_phased:
                raise ValueError(
                    f"Sample node {sid} is not in global_phased_node_ids"
                )
        return sample_sets

    def diversity(
        self,
        sample_sets=None,
        *,
        windows=None,
        mode="site",
        span_normalise=True,
    ):
        """
        Compute diversity across contigs.

        For ``windows=None``, this computes per-contig diversity and combines
        values by contig span so that results match running ``diversity`` on
        ``to_ts(ta)`` for the same arguments.

        Notes
        -----
        ``windows`` is currently unimplemented and will raise
        :class:`NotImplementedError` if provided.
        """
        if windows is not None:
            raise NotImplementedError(
                "ta.stats.diversity does not yet support windows"
            )

        if self._ta.num_contigs == 0:
            raise ValueError("Cannot compute diversity on an empty TreesAssemblage")

        effective_sample_sets = self._resolve_sample_sets(sample_sets)

        values = []
        spans = []
        for ts in self._ta.values():
            values.append(
                ts.diversity(
                    sample_sets=effective_sample_sets,
                    windows=None,
                    mode=mode,
                    span_normalise=span_normalise,
                )
            )
            spans.append(ts.sequence_length)

        spans = np.asarray(spans, dtype=float)
        values = [np.asarray(v, dtype=float) for v in values]

        if span_normalise:
            # Weighted average of per-span-normalised values.
            numer = np.zeros_like(values[0], dtype=float)
            for v, s in zip(values, spans):
                numer = numer + v * s
            out = numer / spans.sum()
        else:
            # Unnormalised values are additive across contigs.
            out = np.zeros_like(values[0], dtype=float)
            for v in values:
                out = out + v

        if np.ndim(out) == 0:
            return float(out)
        return out


def _iter_sample_ids(sample_sets):
    """Yield sample IDs from list[int] or list[list[int]] sample_sets values."""
    for sample_set in sample_sets:
        if np.isscalar(sample_set):
            yield int(sample_set)
        else:
            for sid in sample_set:
                yield int(sid)
