"""
Core TreesAssemblage class for tskit_multichrom.
"""

import collections

import tskit

from .flags import CONTIG_METADATA_KEY, NODE_IS_SHARED

# Named tuple representing a contig's key in the assemblage dictionary.
# - index: ordering integer (not required to be consecutive)
# - id: unique integer identifier
# - symbol: string identifier (e.g. "chrX")
# - type: string type (e.g. "A" for autosome, "X", "Y", ...)
ContigKey = collections.namedtuple("ContigKey", ["index", "id", "symbol", "type"])


class TreesAssemblage:
    """
    A collection of tree sequences representing multiple contigs (chromosomes).

    Each tree sequence is stored under a :class:`ContigKey` namedtuple.
    The assemblage enforces a set of consistency requirements across the
    constituent tree sequences (see :meth:`_validate`).

    Access contigs via :meth:`contig` (by id, symbol, or index).

    Parameters
    ----------
    tree_sequences : dict[ContigKey, tskit.TreeSequence]
        Mapping of contig keys to tree sequences.
    skip_validation : bool
        If True, skip validation on construction. Defaults to False.
    """

    def __init__(self, tree_sequences, *, skip_validation=False):
        if not isinstance(tree_sequences, dict):
            raise TypeError("tree_sequences must be a dict")
        for key, ts in tree_sequences.items():
            if not isinstance(key, ContigKey):
                raise TypeError(f"Keys must be ContigKey instances, got {type(key)}")
            if not isinstance(ts, tskit.TreeSequence):
                raise TypeError(
                    f"Values must be tskit.TreeSequence instances, got {type(ts)}"
                )
        self._tree_sequences = dict(tree_sequences)
        if not skip_validation:
            self._validate()
        self._build_cache()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate(self):
        """Validate consistency requirements across all tree sequences."""
        if len(self._tree_sequences) == 0:
            return

        keys = list(self._tree_sequences.keys())
        tss = list(self._tree_sequences.values())

        # Keys: index and id must be unique; (index, id, symbol) must be unique
        indexes = [k.index for k in keys]
        ids = [k.id for k in keys]
        symbols = [k.symbol for k in keys]

        if len(set(indexes)) != len(indexes):
            raise ValueError("Contig 'index' values must be unique within an assemblage")
        if len(set(ids)) != len(ids):
            raise ValueError("Contig 'id' values must be unique within an assemblage")
        if len(set(symbols)) != len(symbols):
            raise ValueError(
                "Contig 'symbol' values must be unique within an assemblage"
            )
        for k in keys:
            if not isinstance(k.index, int) or k.index < 0:
                raise ValueError(
                    f"Contig index must be a non-negative integer, got {k.index!r}"
                )
            if not isinstance(k.id, int) or k.id < 0:
                raise ValueError(
                    f"Contig id must be a non-negative integer, got {k.id!r}"
                )

        # Check ContigKey matches top-level metadata of each tree sequence
        for key, ts in self._tree_sequences.items():
            self._check_contig_metadata(key, ts)

        # Reference tree sequence (first in index order)
        ref_ts = tss[0]

        # Individual tables must be identical
        for ts in tss[1:]:
            if not ref_ts.tables.individuals.equals(ts.tables.individuals):
                raise ValueError(
                    "Individual tables must be identical across all contigs in an "
                    "assemblage"
                )

        # Population tables must be identical
        for ts in tss[1:]:
            if not ref_ts.tables.populations.equals(ts.tables.populations):
                raise ValueError(
                    "Population tables must be identical across all contigs in an "
                    "assemblage"
                )

        # Migration tables must be empty
        for key, ts in self._tree_sequences.items():
            if ts.num_migrations > 0:
                raise ValueError(
                    f"Migration tables must be empty; contig {key.symbol!r} has "
                    f"{ts.num_migrations} migrations"
                )

        # time_units must be identical
        ref_time_units = ref_ts.time_units
        for ts in tss[1:]:
            if ts.time_units != ref_time_units:
                raise ValueError(
                    "time_units must be identical across all contigs in an assemblage"
                )

        # Node metadata schemas must be identical
        ref_node_schema = ref_ts.tables.nodes.metadata_schema
        for key, ts in self._tree_sequences.items():
            if ts.tables.nodes.metadata_schema != ref_node_schema:
                raise ValueError(
                    f"Node metadata schemas must be identical across all contigs; "
                    f"contig {key.symbol!r} has a different schema"
                )

        # Check IS_SHARED node identity: any node with IS_SHARED set must be
        # identical (flags, time, individual, population, metadata) across all
        # tree sequences that contain that node ID.
        self._validate_shared_nodes()

    @staticmethod
    def _check_contig_metadata(key, ts):
        """Check that the ContigKey matches the contig metadata in the tree sequence."""
        meta = ts.metadata
        if not isinstance(meta, dict) or CONTIG_METADATA_KEY not in meta:
            raise ValueError(
                f"Tree sequence for contig {key.symbol!r} must have a "
                f"'{CONTIG_METADATA_KEY}' key in its top-level metadata"
            )
        contig_meta = meta[CONTIG_METADATA_KEY]
        for field in ("index", "id", "symbol", "type"):
            if field not in contig_meta:
                raise ValueError(
                    f"Contig metadata for {key.symbol!r} is missing field {field!r}"
                )
        if contig_meta["index"] != key.index:
            raise ValueError(
                f"ContigKey.index ({key.index}) does not match metadata 'index' "
                f"({contig_meta['index']}) for contig {key.symbol!r}"
            )
        if contig_meta["id"] != key.id:
            raise ValueError(
                f"ContigKey.id ({key.id}) does not match metadata 'id' "
                f"({contig_meta['id']}) for contig {key.symbol!r}"
            )
        if contig_meta["symbol"] != key.symbol:
            raise ValueError(
                f"ContigKey.symbol ({key.symbol!r}) does not match metadata "
                f"'symbol' ({contig_meta['symbol']!r})"
            )
        if contig_meta["type"] != key.type:
            raise ValueError(
                f"ContigKey.type ({key.type!r}) does not match metadata "
                f"'type' ({contig_meta['type']!r})"
            )

    def _validate_shared_nodes(self):
        """
        Validate IS_SHARED node identity requirements.

        Any node with IS_SHARED flag set must be identical (same flags, time,
        individual, population, metadata) across all tree sequences that contain
        that node ID.
        """
        # Collect all shared node IDs and their representative rows
        shared_nodes = {}  # node_id -> (representative_row, contig_symbol)

        for key, ts in self._tree_sequences.items():
            nodes = ts.tables.nodes
            for node_id in range(ts.num_nodes):
                row = nodes[node_id]
                if row.flags & NODE_IS_SHARED:
                    if node_id in shared_nodes:
                        ref_row, ref_symbol = shared_nodes[node_id]
                        # Must be identical
                        if not _node_rows_equal(ref_row, row):
                            raise ValueError(
                                f"Node {node_id} has IS_SHARED flag set but differs "
                                f"between contig {ref_symbol!r} and {key.symbol!r}: "
                                f"{ref_row} vs {row}"
                            )
                    else:
                        shared_nodes[node_id] = (row, key.symbol)

    # ------------------------------------------------------------------
    # Caching
    # ------------------------------------------------------------------

    def _build_cache(self):
        """Build the internal cache for the assemblage."""
        self._cache = {}

        # Sort keys by index
        sorted_keys = sorted(self._tree_sequences.keys(), key=lambda k: k.index)
        self._sorted_keys = sorted_keys

        # Total sequence length (sum of all contig lengths)
        total_length = sum(
            ts.sequence_length for ts in self._tree_sequences.values()
        )
        self._cache["total_sequence_length"] = total_length

        # Cross-phased node IDs: nodes with IS_SHARED flag in ALL tree sequences
        # that contain that node ID
        self._cache["cross_phased_node_ids"] = self._compute_cross_phased_nodes()

        # Count nonglobal sample nodes: sample nodes that are NOT cross-phased
        # in at least one contig
        all_sample_ids = set()
        for ts in self._tree_sequences.values():
            for s in ts.samples():
                all_sample_ids.add(s)
        cross_phased = self._cache["cross_phased_node_ids"]
        nonglobal = all_sample_ids - cross_phased
        self._cache["nonglobal_sample_node_count"] = len(nonglobal)
        self._cache["is_partial_sample_arg"] = len(nonglobal) > 0

        # Index by id and symbol for fast lookup
        self._by_id = {k.id: k for k in self._tree_sequences}
        self._by_symbol = {k.symbol: k for k in self._tree_sequences}

    def _compute_cross_phased_nodes(self):
        """
        Return a set of node IDs that have IS_SHARED flag set in every tree
        sequence that contains them.

        A node ID is considered "cross-phased" if:
        - It exists in at least one tree sequence, AND
        - In every tree sequence that contains it, it has the IS_SHARED flag set.
        """
        if not self._tree_sequences:
            return set()

        # node_id -> count of TS where it appears with IS_SHARED set
        shared_count = {}
        # node_id -> count of TS where it appears (at all)
        total_count = {}

        for ts in self._tree_sequences.values():
            for node_id in range(ts.num_nodes):
                total_count[node_id] = total_count.get(node_id, 0) + 1
                if ts.tables.nodes[node_id].flags & NODE_IS_SHARED:
                    shared_count[node_id] = shared_count.get(node_id, 0) + 1

        # A node is cross-phased if shared_count == total_count for that node
        return {
            nid
            for nid, cnt in total_count.items()
            if shared_count.get(nid, 0) == cnt
        }

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def tree_sequences(self):
        """dict[ContigKey, tskit.TreeSequence]: The underlying tree sequences."""
        return dict(self._tree_sequences)

    @property
    def contigs(self):
        """list[ContigKey]: Contig keys sorted by index."""
        return list(self._sorted_keys)

    @property
    def total_sequence_length(self):
        """float: Sum of sequence_length of all contigs."""
        return self._cache["total_sequence_length"]

    @property
    def cross_phased_node_ids(self):
        """set[int]: Node IDs that are shared (cross-phased) in all contigs."""
        return set(self._cache["cross_phased_node_ids"])

    @property
    def nonglobal_sample_node_count(self):
        """int: Number of sample nodes not shared in all contigs."""
        return self._cache["nonglobal_sample_node_count"]

    @property
    def is_partial_sample_arg(self):
        """bool: True if there are any nonglobal sample nodes."""
        return self._cache["is_partial_sample_arg"]

    @property
    def num_contigs(self):
        """int: Number of contigs in the assemblage."""
        return len(self._tree_sequences)

    # ------------------------------------------------------------------
    # Contig access
    # ------------------------------------------------------------------

    def contig(self, id_or_symbol):
        """
        Return the tree sequence for a contig by its symbol (str), integer id,
        or integer index.

        If a string is provided, it is treated as the contig symbol (e.g.
        ``"chrX"``). If an integer is provided it is first matched against the
        contig ``id`` field, and then against the contig ``index`` field.

        Parameters
        ----------
        id_or_symbol : str or int

        Returns
        -------
        tskit.TreeSequence
        """
        if isinstance(id_or_symbol, str):
            if id_or_symbol not in self._by_symbol:
                raise KeyError(
                    f"No contig with symbol {id_or_symbol!r}. "
                    f"Available: {list(self._by_symbol)}"
                )
            return self._tree_sequences[self._by_symbol[id_or_symbol]]

        # Integer lookup: try id first, then index
        if id_or_symbol in self._by_id:
            return self._tree_sequences[self._by_id[id_or_symbol]]
        for key in self._tree_sequences:
            if key.index == id_or_symbol:
                return self._tree_sequences[key]
        raise KeyError(
            f"No contig with id or index {id_or_symbol!r}. "
            f"Available ids: {list(self._by_id)}"
        )

    def __len__(self):
        return len(self._tree_sequences)

    def __iter__(self):
        """Iterate over ContigKeys in index order."""
        return iter(self._sorted_keys)

    def __getitem__(self, key):
        """Return the tree sequence for a ContigKey."""
        return self._tree_sequences[key]

    def __repr__(self):
        symbols = [k.symbol for k in self._sorted_keys]
        return (
            f"TreesAssemblage(contigs={symbols!r}, "
            f"total_length={self.total_sequence_length})"
        )

    # ------------------------------------------------------------------
    # Subsetting
    # ------------------------------------------------------------------

    def subset(self, *, symbols=None, type=None, types=None, ids=None, indexes=None):
        """
        Create a new :class:`TreesAssemblage` from a subset of contigs.

        Subsetting is cheap: it does not copy the underlying tree sequences.
        Only the cache is recomputed.

        Parameters
        ----------
        symbols : list[str], optional
            Keep only contigs with these symbols.
        type : str, optional
            Keep only contigs whose type matches this single string value.
        types : list[str], optional
            Keep only contigs with these types (e.g. ``["A"]`` for autosomes).
            If both ``type`` and ``types`` are given, they are combined.
        ids : list[int], optional
            Keep only contigs with these id values.
        indexes : list[int], optional
            Keep only contigs with these index values.

        Returns
        -------
        TreesAssemblage
        """
        keep = set(self._tree_sequences.keys())

        if symbols is not None:
            symbols_set = set(symbols)
            keep &= {k for k in keep if k.symbol in symbols_set}

        # Merge `type` (single) and `types` (list) into one filter
        types_filter = None
        if type is not None:
            types_filter = {type}
        if types is not None:
            types_filter = (types_filter or set()) | set(types)
        if types_filter is not None:
            keep &= {k for k in keep if k.type in types_filter}

        if ids is not None:
            ids_set = set(ids)
            keep &= {k for k in keep if k.id in ids_set}
        if indexes is not None:
            indexes_set = set(indexes)
            keep &= {k for k in keep if k.index in indexes_set}

        new_ts = {k: self._tree_sequences[k] for k in keep}
        return TreesAssemblage(new_ts, skip_validation=True)

    # ------------------------------------------------------------------
    # Reindexing
    # ------------------------------------------------------------------

    def reindex(self, order=None):
        """
        Return a new :class:`TreesAssemblage` with contigs reindexed from 0..N-1.

        This makes copies of the underlying tree sequences with updated metadata.

        Parameters
        ----------
        order : list[str or int], optional
            List of symbols or ids specifying the desired order.
            If None, uses current index order.

        Returns
        -------
        TreesAssemblage
        """
        if order is None:
            ordered_keys = self._sorted_keys
        else:
            ordered_keys = []
            for item in order:
                if isinstance(item, str):
                    if item not in self._by_symbol:
                        raise KeyError(f"No contig with symbol {item!r}")
                    ordered_keys.append(self._by_symbol[item])
                elif isinstance(item, int):
                    if item not in self._by_id:
                        raise KeyError(f"No contig with id {item!r}")
                    ordered_keys.append(self._by_id[item])
                else:
                    raise TypeError(
                        f"order items must be str (symbol) or int (id), got {type(item)}"
                    )
            if len(set(ordered_keys)) != len(self._tree_sequences):
                raise ValueError(
                    "order must specify all contigs exactly once"
                )

        new_ts = {}
        for new_index, key in enumerate(ordered_keys):
            ts = self._tree_sequences[key]
            tables = ts.dump_tables()
            meta = dict(ts.metadata)
            contig_meta = dict(meta[CONTIG_METADATA_KEY])
            contig_meta["index"] = new_index
            meta[CONTIG_METADATA_KEY] = contig_meta
            tables.metadata = meta
            new_key = ContigKey(
                index=new_index,
                id=key.id,
                symbol=key.symbol,
                type=key.type,
            )
            new_ts[new_key] = tables.tree_sequence()

        return TreesAssemblage(new_ts)


# Backwards-compatible alias
TreesArchive = TreesAssemblage


# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------


def _node_rows_equal(row1, row2):
    """Compare two node table rows for equality (flags, time, individual, pop, meta)."""
    return (
        row1.flags == row2.flags
        and row1.time == row2.time
        and row1.individual == row2.individual
        and row1.population == row2.population
        and row1.metadata == row2.metadata
    )


def make_contig_key(contig_meta):
    """
    Create a :class:`ContigKey` from a contig metadata dict.

    Parameters
    ----------
    contig_meta : dict
        Dictionary with keys 'index', 'id', 'symbol', 'type'.

    Returns
    -------
    ContigKey
    """
    return ContigKey(
        index=int(contig_meta["index"]),
        id=int(contig_meta["id"]),
        symbol=str(contig_meta["symbol"]),
        type=str(contig_meta["type"]),
    )


def make_permissive_contig_schema():
    """
    Return a permissive :class:`tskit.MetadataSchema` suitable for contig
    top-level metadata.

    Uses :meth:`tskit.MetadataSchema.permissive_json` as the base schema so
    that any additional top-level keys are tolerated. The presence of the
    required ``'contig'`` key is validated separately when a
    :class:`TreesAssemblage` is constructed (see
    :meth:`TreesAssemblage._check_contig_metadata`).
    """
    return tskit.MetadataSchema.permissive_json()
