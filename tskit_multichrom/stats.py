"""Statistics helpers for TreeSequenceDictionary."""

import numpy as np


class TreeSequenceDictionaryStats:
    """Statistics namespace for :class:`tskit_multichrom.core.TreeSequenceDictionary`."""

    def __init__(self, tsd):
        self._tsd = tsd

    def _resolve_sample_sets(self, sample_sets):
        """
        Resolve and validate ``sample_sets`` for multi-contig statistics.

        Rules:
        - If ``sample_sets is None`` and there are no nonglobal samples,
          default to one sample set containing all ``global_phased_node_ids``.
        - If ``sample_sets is None`` and nonglobal samples exist, raise.
        - If ``sample_sets`` is provided, every sample ID must be globally phased.
        """
        global_phased = set(self._tsd.global_phased_node_ids)

        if sample_sets is None:
            if self._tsd.nonglobal_sample_node_count != 0:
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
        ``to_ts(tsd)`` for the same arguments.

        Notes
        -----
        ``windows`` is currently unimplemented and will raise
        :class:`NotImplementedError` if provided.
        """
        if windows is not None:
            raise NotImplementedError(
                "tsd.stats.diversity does not yet support windows"
            )

        if self._tsd.num_contigs == 0:
            raise ValueError("Cannot compute diversity on an empty TreeSequenceDictionary")

        effective_sample_sets = self._resolve_sample_sets(sample_sets)

        values = []
        spans = []
        for ts in self._tsd.values():
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

    def pca(
        self,
        num_components,
        *,
        windows=None,
        samples=None,
        individuals=None,
        time_windows=None,
        mode="branch",
        centre=True,
        num_iterations=5,
        num_oversamples=None,
        random_seed=None,
        range_sketch=None,
    ):
        """
        Perform PCA across all contigs.

        Internally calls ``tsd.to_ts()`` to obtain a single combined tree
        sequence and then delegates to :meth:`tskit.TreeSequence.pca`.

        When ``samples`` is provided (or defaulted) every sample node ID must
        be globally phased, consistent with the requirements of
        :meth:`diversity`. When ``individuals`` is provided, all contigs must
        share the same chromosome type so that per-individual ploidy is
        consistent across contigs.

        Parameters
        ----------
        num_components : int
            Number of principal components to return.
        windows : list, optional
            Genomic windows; coordinates refer to the combined tree sequence
            returned by ``tsd.to_ts()``.
        samples : array_like, optional
            Sample node IDs. Must be globally phased. Mutually exclusive with
            ``individuals``.
        individuals : array_like, optional
            Individual IDs forwarded to ``ts.pca()``. All contigs must share
            the same chromosome type. Mutually exclusive with ``samples``.
        time_windows, mode, centre, num_iterations, num_oversamples,
        random_seed, range_sketch :
            Forwarded verbatim to :meth:`tskit.TreeSequence.pca`.

        Returns
        -------
        tskit.PCAResult

        Notes
        -----
        TODO: A more efficient implementation could combine per-contig PCA
        results mathematically (via a block-structured GRM), avoiding the cost
        of materialising ``to_ts()`` for large assemblages.
        """
        if samples is not None and individuals is not None:
            raise ValueError("Cannot specify both samples and individuals")

        ts = self._tsd.to_ts()

        if individuals is not None:
            types = {key.type for key in self._tsd.contigs}
            if len(types) > 1:
                raise ValueError(
                    "pca with individuals requires all contigs to be the same"
                    f" type; found types: {sorted(types)}"
                )
            # tskit's pca requires sample node IDs to be contiguous [0, n).
            # to_ts() may interleave ancestral nodes among non-shared sample
            # nodes, so simplify first to compact sample IDs.
            ts = ts.simplify(samples=ts.samples(), record_provenance=False)
            return ts.pca(
                num_components=num_components,
                windows=windows,
                individuals=individuals,
                time_windows=time_windows,
                mode=mode,
                centre=centre,
                num_iterations=num_iterations,
                num_oversamples=num_oversamples,
                random_seed=random_seed,
                range_sketch=range_sketch,
            )

        if samples is None:
            if self._tsd.is_nonglobal_sample_arg:
                raise ValueError(
                    "individuals or samples must be provided when nonglobal"
                    " sample nodes are present"
                )
            effective_samples = np.asarray(
                sorted(self._tsd.global_phased_node_ids), dtype=np.int32
            )
        else:
            global_phased = set(self._tsd.global_phased_node_ids)
            for sid in samples:
                if int(sid) not in global_phased:
                    raise ValueError(
                        f"Sample node {sid} is not in global_phased_node_ids"
                    )
            effective_samples = samples

        return ts.pca(
            num_components=num_components,
            windows=windows,
            samples=effective_samples,
            time_windows=time_windows,
            mode=mode,
            centre=centre,
            num_iterations=num_iterations,
            num_oversamples=num_oversamples,
            random_seed=random_seed,
            range_sketch=range_sketch,
        )


def _iter_sample_ids(sample_sets):
    """Yield sample IDs from list[int] or list[list[int]] sample_sets values."""
    for sample_set in sample_sets:
        if np.isscalar(sample_set):
            yield int(sample_set)
        else:
            for sid in sample_set:
                yield int(sid)

