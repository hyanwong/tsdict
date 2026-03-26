"""
Conversion functions between a :class:`TreeSequenceDictionary` and a single
:class:`tskit.TreeSequence`.

Overview
--------
**to_ts** merges all contigs in a :class:`TreeSequenceDictionary` into a single
:class:`tskit.TreeSequence` by concatenating their genomic coordinate systems,
placing contigs end-to-end.
Shared nodes (IS_SHARED flag) retain their IDs; non-shared nodes get new IDs.

**from_ts** reverses this process, splitting the single :class:`tskit.TreeSequence`
back into a :class:`TreeSequenceDictionary` using the top-level metadata array that
records per-contig ``sequence_length``, ``num_nodes``, and ``contig``.

**from_slim** converts a set of SLiM-style tree sequences (where all nodes
are implicitly shared) into a :class:`TreeSequenceDictionary`.
"""

import copy
import warnings
import json

import numpy as np
import tskit

from ._version import tsdict_version
from .core import ContigKey, TreeSequenceDictionary, make_contig_key
from .flags import CONTIG_METADATA_KEY, NODE_IS_SHARED

# Top-level metadata key used in the single-TS representation
_ARCHIVE_META_KEY = "contigs"


def to_ts(tsd, record_provenance=True):
    """
    Merge a :class:`TreeSequenceDictionary` into a single :class:`tskit.TreeSequence`.

    Contigs are placed end-to-end in order of their ``index`` values.
    Shared nodes (those with :data:`~tsdict.flags.NODE_IS_SHARED` set)
    retain their original node IDs in the merged tree sequence. Non-shared nodes
    receive new IDs, assigned after all shared nodes have been placed.

    Site and mutation metadata schemas must be identical across all contigs,
    otherwise a :class:`ValueError` is raised.

    Parameters
    ----------
    tsd : TreeSequenceDictionary
    record_provenance : bool
        Whether to append a provenance record for this ``to_ts`` call.

    Returns
    -------
    tskit.TreeSequence
    """
    if tsd.num_contigs == 0:
        raise ValueError("Cannot convert an empty archive to a tree sequence")

    sorted_keys = tsd.contigs  # sorted by index

    # ------------------------------------------------------------------
    # 1. Validate that site and mutation schemas are identical
    # ------------------------------------------------------------------
    ref_site_schema = tsd[sorted_keys[0]].tables.sites.metadata_schema
    ref_mut_schema = tsd[sorted_keys[0]].tables.mutations.metadata_schema
    for key in sorted_keys[1:]:
        ts = tsd[key]
        if ts.tables.sites.metadata_schema != ref_site_schema:
            raise ValueError(
                f"Site metadata schemas differ between contigs: cannot merge. "
                f"Contig {key.symbol!r} has a different schema."
            )
        if ts.tables.mutations.metadata_schema != ref_mut_schema:
            raise ValueError(
                f"Mutation metadata schemas differ between contigs: cannot merge. "
                f"Contig {key.symbol!r} has a different schema."
            )

    # ------------------------------------------------------------------
    # 2. Build the shared_node_array
    #    Each entry: (node_id, first_chrom_index, vacant_bitflags)
    #    vacant_bitflags: bit i set if the contig with metadata 'index' i is absent
    # ------------------------------------------------------------------
    # Node IDs with IS_SHARED set in at least one contig.
    shared_node_ids = tsd.shared_node_ids

    # For each shared node ID, build vacant_bitflags and find first non-vacant contig
    shared_node_array = []  # list of (node_id, first_contig_key, vacant_bitflags)
    for nid in sorted(shared_node_ids):
        # vacant_bitflags: bit = contig 'index' value
        vacant_bitflags = 0
        first_key = None
        for key in sorted_keys:
            ts = tsd[key]
            if nid < ts.num_nodes and (ts.tables.nodes[nid].flags & NODE_IS_SHARED):
                if first_key is None:
                    first_key = key
            else:
                # This contig is vacant for this node
                vacant_bitflags |= 1 << key.index
        if first_key is None:
            continue  # shouldn't happen
        shared_node_array.append((nid, first_key, vacant_bitflags))

    # ------------------------------------------------------------------
    # 3. Build new TableCollection
    # ------------------------------------------------------------------
    ref_ts = tsd[sorted_keys[0]]
    new_tc = tskit.TableCollection(sequence_length=tsd.total_sequence_length)

    # Copy individual, population tables from reference
    new_tc.individuals.replace_with(ref_ts.tables.individuals)
    new_tc.populations.replace_with(ref_ts.tables.populations)

    # Always copy provenance from reference (first in index order).
    new_tc.provenances.replace_with(ref_ts.tables.provenances)

    # Set metadata schemas
    new_tc.nodes.metadata_schema = ref_ts.tables.nodes.metadata_schema
    new_tc.sites.metadata_schema = ref_site_schema
    new_tc.mutations.metadata_schema = ref_mut_schema

    # ------------------------------------------------------------------
    # 4. Build top-level metadata for the merged TS
    # ------------------------------------------------------------------
    contigs_meta = []
    for key in sorted_keys:
        ts = tsd[key]
        entry = {
            "sequence_length": ts.sequence_length,
            "num_nodes": ts.num_nodes,
            CONTIG_METADATA_KEY: dict(ts.metadata.get(CONTIG_METADATA_KEY, {})),
        }
        contigs_meta.append(entry)

    permissive_schema = {
        "codec": "json",
        "type": "object",
        "additionalProperties": True,
    }
    new_tc.metadata_schema = tskit.MetadataSchema(permissive_schema)
    new_tc.metadata = {_ARCHIVE_META_KEY: contigs_meta}

    if record_provenance:
        provenance_record = {
            "software": {
                "name": "tsdict",
                "version": tsdict_version,
            },
            "parameters": {
                "command": "to_ts",
                "num_contigs": len(sorted_keys),
            },
        }
        new_tc.provenances.add_row(record=json.dumps(provenance_record))

    # ------------------------------------------------------------------
    # 5. Iterate through contigs and build node/edge/site/mutation tables
    # ------------------------------------------------------------------
    start_position = 0.0
    shared_node_counter = 0

    # node_map[contig_index][old_node_id] = new_node_id
    all_node_maps = []

    for contig_idx, key in enumerate(sorted_keys):
        ts = tsd[key]
        node_map = np.full(ts.num_nodes, tskit.NULL, dtype=np.int64)

        # Add nodes from this contig, inserting shared nodes at correct positions
        for old_node_id in range(ts.num_nodes):
            # Check if we need to insert a shared node here
            while (
                shared_node_counter < len(shared_node_array)
                and shared_node_array[shared_node_counter][0]
                == len(new_tc.nodes)
            ):
                shared_nid, first_key, vacant_bitflags = shared_node_array[
                    shared_node_counter
                ]
                # Append node from the first non-vacant contig
                src_ts = tsd[first_key]
                src_row = src_ts.tables.nodes[shared_nid]

                # Optionally store vacant_bitflags in metadata
                node_meta = _try_add_vacant_flags(
                    src_row.metadata,
                    vacant_bitflags,
                    src_ts.tables.nodes.metadata_schema,
                )
                new_tc.nodes.append(
                    tskit.NodeTableRow(
                        flags=src_row.flags,
                        time=src_row.time,
                        population=src_row.population,
                        individual=src_row.individual,
                        metadata=node_meta,
                    )
                )
                shared_node_counter += 1

            # Now handle old_node_id from this contig
            row = ts.tables.nodes[old_node_id]
            if (row.flags & NODE_IS_SHARED) and old_node_id in shared_node_ids:
                # This node will appear in new_tc with the same ID
                # Its position in new_tc is exactly old_node_id
                # (assuming shared node IDs are low and consecutive enough)
                # We need to find where it actually ended up in new_tc
                # Since we insert shared nodes in order, we can look it up
                # by scanning for the node in new_tc — but more efficiently,
                # shared nodes with IS_SHARED should map to their own ID
                node_map[old_node_id] = old_node_id
            else:
                # Non-shared node: gets a new ID
                new_id = len(new_tc.nodes)
                node_map[old_node_id] = new_id
                new_tc.nodes.append(
                    tskit.NodeTableRow(
                        flags=row.flags,
                        time=row.time,
                        population=row.population,
                        individual=row.individual,
                        metadata=row.metadata,
                    )
                )

        # Flush any remaining shared nodes that fall after all nodes in this contig
        # (this handles shared nodes that are beyond num_nodes for this contig)

        all_node_maps.append(node_map)

        # Add edges: shift coordinates and remap node IDs
        edges = ts.tables.edges
        if len(edges) > 0:
            new_tc.edges.append_columns(
                left=edges.left + start_position,
                right=edges.right + start_position,
                parent=node_map[edges.parent].astype(np.int32),
                child=node_map[edges.child].astype(np.int32),
                metadata=edges.metadata,
                metadata_offset=edges.metadata_offset,
            )

        # Add sites: shift positions
        sites = ts.tables.sites
        site_id_offset = len(new_tc.sites)
        if len(sites) > 0:
            new_tc.sites.append_columns(
                position=sites.position + start_position,
                ancestral_state=sites.ancestral_state,
                ancestral_state_offset=sites.ancestral_state_offset,
                metadata=sites.metadata,
                metadata_offset=sites.metadata_offset,
            )

        # Add mutations: remap node IDs and site IDs, fix parent references
        mutations = ts.tables.mutations
        mut_id_offset_before = len(new_tc.mutations)
        if len(mutations) > 0:
            new_parents = np.where(
                mutations.parent == tskit.NULL,
                tskit.NULL,
                mutations.parent + mut_id_offset_before,
            ).astype(np.int32)
            new_tc.mutations.append_columns(
                site=mutations.site + site_id_offset,
                node=node_map[mutations.node].astype(np.int32),
                time=mutations.time,
                derived_state=mutations.derived_state,
                derived_state_offset=mutations.derived_state_offset,
                parent=new_parents,
                metadata=mutations.metadata,
                metadata_offset=mutations.metadata_offset,
            )

        start_position += ts.sequence_length

    # Flush any remaining shared nodes (those with IDs beyond all contig node counts)
    while shared_node_counter < len(shared_node_array):
        shared_nid, first_key, vacant_bitflags = shared_node_array[shared_node_counter]
        # Pad up to this node ID if needed
        while len(new_tc.nodes) < shared_nid:
            # Add a placeholder node (should not normally happen)
            warnings.warn(
                f"Padding node table to reach shared node ID {shared_nid}",
                stacklevel=2,
            )
            new_tc.nodes.append(
                tskit.NodeTableRow(flags=0, time=0.0, population=tskit.NULL, individual=tskit.NULL)
            )
        src_ts = tsd[first_key]
        src_row = src_ts.tables.nodes[shared_nid]
        node_meta = _try_add_vacant_flags(
            src_row.metadata,
            vacant_bitflags,
            src_ts.tables.nodes.metadata_schema,
        )
        new_tc.nodes.append(
            tskit.NodeTableRow(
                flags=src_row.flags,
                time=src_row.time,
                population=src_row.population,
                individual=src_row.individual,
                metadata=node_meta,
            )
        )
        shared_node_counter += 1

    new_tc.sort()
    return new_tc.tree_sequence()


def from_ts(ts, record_provenance=True):
    """
    Convert a single merged :class:`tskit.TreeSequence` back into a
    :class:`TreeSequenceDictionary`.

    The input tree sequence must have been created by :func:`to_ts`
    (or equivalent), with top-level metadata containing a
    ``'tsdict_contigs'`` array.

    Parameters
    ----------
    ts : tskit.TreeSequence
    record_provenance : bool
        Whether to record provenance while splitting the merged tree sequence.

    Returns
    -------
    TreeSequenceDictionary
    """
    meta = ts.metadata
    if not isinstance(meta, dict) or _ARCHIVE_META_KEY not in meta:
        raise ValueError(
            f"Tree sequence does not have '{_ARCHIVE_META_KEY}' in top-level metadata. "
            "Was it created by tsdict.to_ts()?"
        )

    contigs_meta = meta[_ARCHIVE_META_KEY]

    # Validate cumulative lengths
    sequence_lengths = [c["sequence_length"] for c in contigs_meta]
    cumsum = [0.0] + list(np.cumsum(sequence_lengths))

    if abs(cumsum[-1] - ts.sequence_length) > 1e-10:
        raise ValueError(
            f"Sum of contig sequence lengths ({cumsum[-1]}) does not match "
            f"tree sequence length ({ts.sequence_length})"
        )

    # Check that index values are increasing in order
    indexes = [c[CONTIG_METADATA_KEY]["index"] for c in contigs_meta]
    if indexes != sorted(indexes):
        raise ValueError(
            f"Contig indexes in metadata must be in increasing order; got {indexes}"
        )

    # Build a boolean array of nodes used so far
    node_used = np.zeros(ts.num_nodes, dtype=bool)

    # Determine which nodes have IS_SHARED set
    is_shared_global = (ts.tables.nodes.flags & NODE_IS_SHARED).astype(bool)

    result = {}

    for i, contig_entry in enumerate(contigs_meta):
        contig_meta_dict = contig_entry[CONTIG_METADATA_KEY]
        seq_len = contig_entry["sequence_length"]
        n_nodes = contig_entry["num_nodes"]
        key = make_contig_key(contig_meta_dict)

        left = cumsum[i]
        right = cumsum[i + 1]

        # Determine which nodes to keep for this contig:
        # - non-shared nodes that have been used so far (sequentially assigned)
        # - shared nodes that are non-vacant for this contig
        is_shared_for_contig = _get_shared_for_contig(
            ts, is_shared_global, contig_meta_dict.get("index", i)
        )
        non_shared_used = node_used & ~is_shared_global
        candidate = non_shared_used | is_shared_for_contig
        # Take the first n_nodes True values as node IDs
        candidate_ids = np.where(candidate)[0]
        if len(candidate_ids) < n_nodes:
            # Need more nodes — take next available non-shared nodes
            all_ids = np.arange(ts.num_nodes)
            not_yet_assigned = ~node_used & ~is_shared_global
            extra = all_ids[not_yet_assigned][: n_nodes - len(candidate_ids)]
            candidate_ids = np.sort(np.concatenate([candidate_ids, extra]))
        else:
            candidate_ids = candidate_ids[:n_nodes]

        node_used[candidate_ids] = True

        # Extract contig tree sequence using keep_intervals then shift
        contig_ts = (
            ts.keep_intervals(
                [[left, right]],
                simplify=False,
                record_provenance=False,
            )
            .shift(
                -left,
                sequence_length=seq_len,
                record_provenance=False,
            )
        )

        # Subset to the correct nodes
        # We need to keep candidate_ids in the new (shifted) tree sequence
        contig_tc = contig_ts.dump_tables()
        contig_tc.subset(candidate_ids, record_provenance=False)

        # Set per-contig metadata schema and metadata
        contig_tc.metadata_schema = tskit.MetadataSchema.permissive_json()
        contig_tc.metadata = {CONTIG_METADATA_KEY: contig_meta_dict}

        if record_provenance:
            provenance_record = {
                "software": {
                    "name": "tsdict",
                    "version": tsdict_version,
                },
                "parameters": {
                    "command": "from_ts",
                    "contig_index": int(contig_meta_dict.get("index", i)),
                    "contig_symbol": str(contig_meta_dict.get("symbol", "")),
                },
            }
            contig_tc.provenances.add_row(record=json.dumps(provenance_record))

        result[key] = contig_tc.tree_sequence()

    return TreeSequenceDictionary(result)


def from_slim(tree_sequences, *, slim_metadata_key="SLiM"):
    """
    Convert SLiM-style tree sequences into a :class:`TreeSequenceDictionary`.

    In SLiM tree sequences all nodes are shared across chromosomes, so
    :data:`~tsdict.flags.NODE_IS_SHARED` is set on every node.
    The contig information is read from the SLiM top-level metadata
    ``['SLiM']['this_chromosome']`` and duplicated into the ``'contig'`` key.

    Parameters
    ----------
    tree_sequences : list[tskit.TreeSequence]
        List of SLiM-generated tree sequences, one per chromosome.
    slim_metadata_key : str
        Top-level metadata key for SLiM metadata (default ``'SLiM'``).

    Returns
    -------
    TreeSequenceDictionary
    """
    result = {}
    for ts in tree_sequences:
        meta = ts.metadata
        # When there's no metadata schema, metadata is raw bytes
        if isinstance(meta, (bytes, bytearray)):
            raise ValueError(
                f"Tree sequence does not have '{slim_metadata_key}' in top-level "
                f"metadata; is this a SLiM tree sequence? (metadata is raw bytes, "
                f"no JSON schema found)"
            )
        if slim_metadata_key not in meta:
            raise ValueError(
                f"Tree sequence does not have '{slim_metadata_key}' in top-level "
                f"metadata; is this a SLiM tree sequence?"
            )
        slim_meta = meta[slim_metadata_key]
        if "this_chromosome" not in slim_meta:
            raise ValueError(
                f"SLiM metadata does not have 'this_chromosome'; "
                f"got keys: {list(slim_meta)}"
            )
        chrom = slim_meta["this_chromosome"]
        contig_dict = {
            "index": int(chrom["index"]),
            "id": int(chrom["id"]),
            "symbol": str(chrom["symbol"]),
            "type": str(chrom["type"]),
        }
        key = make_contig_key(contig_dict)

        # Mark all nodes as shared
        tables = ts.dump_tables()
        tables.nodes.flags = tables.nodes.flags | NODE_IS_SHARED

        # Add 'contig' key to top-level metadata
        new_meta = dict(ts.metadata)
        new_meta[CONTIG_METADATA_KEY] = contig_dict

        # Ensure the top-level metadata schema is permissive JSON
        existing_schema = ts.metadata_schema.schema
        if not existing_schema or existing_schema.get("codec") != "json":
            tables.metadata_schema = tskit.MetadataSchema.permissive_json()
        else:
            # Try to add contig key to existing schema (permissive)
            import copy
            schema_dict = copy.deepcopy(existing_schema)
            schema_dict.setdefault("additionalProperties", True)
            if "properties" not in schema_dict:
                schema_dict["properties"] = {}
            schema_dict["properties"][CONTIG_METADATA_KEY] = {
                "type": "object",
                "additionalProperties": True,
            }
            tables.metadata_schema = tskit.MetadataSchema(schema_dict)

        tables.metadata = new_meta
        result[key] = tables.tree_sequence()

    return TreeSequenceDictionary(result)


def from_tree_sequences(
    tree_sequences,
    ids,
    symbols,
    types,
    indexes=None,
    shared_nodes=None,
):
    """
    Create a :class:`TreeSequenceDictionary` from a list of tree sequences without contig metadata.

    This function combines multiple independent tree sequences into a TreeSequenceDictionary,
    automatically annotating them with contig metadata and optionally marking shared nodes.

    Parameters
    ----------
    tree_sequences : list[tskit.TreeSequence]
        Tree sequences to assemble, one per contig.
    ids : list[int]
        Contig IDs (one per tree sequence).
    symbols : list[str]
        Contig symbols, e.g., ``["20", "21"]`` or ``["chrX", "chrY"]``
        (one per tree sequence).
    types : list[str]
        Contig types (one per tree sequence). Examples: ``"A"`` for autosome,
        ``"X"``, ``"Y"``, or ``"MT"``.
    indexes : list[int], optional
        Contig indexes (ordering integers, one per tree sequence). If None,
        defaults to ``[0, 1, 2, ...]``.
    shared_nodes : str or list[int], optional
        Which nodes to mark with :data:`~tsdict.flags.NODE_IS_SHARED`:

        - ``None`` (default): Do not mark any nodes as IS_SHARED.
        - ``"samples"``: Mark all sample nodes as IS_SHARED.
        - list-like of int: Mark the specified node IDs as IS_SHARED across all contigs.

    Returns
    -------
    TreeSequenceDictionary
        A new TreeSequenceDictionary with contigs in index order.

    Raises
    ------
    ValueError
        If the lengths of ``ids``, ``symbols``, ``types``, or ``indexes`` do not match
        the length of ``tree_sequences``, or if ``shared_nodes`` is neither None, the string
        ``"samples"``, nor a list-like of integers.

    Examples
    --------
    >>> tsd = from_tree_sequences(
    ...     [ts20, ts21],
    ...     ids=[20, 21],
    ...     symbols=["20", "21"],
    ...     types=["A", "A"],
    ... )

    Mark all sample nodes as shared:

    >>> tsd = from_tree_sequences(
    ...     [ts20, ts21],
    ...     ids=[20, 21],
    ...     symbols=["20", "21"],
    ...     types=["A", "A"],
    ...     shared_nodes="samples",
    ... )
    """
    if not tree_sequences:
        raise ValueError("tree_sequences list cannot be empty")

    n = len(tree_sequences)
    if len(ids) != n or len(symbols) != n or len(types) != n:
        raise ValueError(
            "ids, symbols, and types must have the same length as tree_sequences"
        )

    if indexes is None:
        indexes = list(range(n))
    elif len(indexes) != n:
        raise ValueError(
            "If provided, indexes must have the same length as tree_sequences"
        )

    # Validate and process shared_nodes parameter
    # mark_strategy will be:
    #   None = don't mark anything
    #   "samples" = mark all sample nodes
    #   set of ints = mark specific node IDs
    if shared_nodes is None:
        mark_strategy = None
    elif isinstance(shared_nodes, str):
        if shared_nodes != "samples":
            raise ValueError(
                f"shared_nodes string must be 'samples', got {shared_nodes!r}"
            )
        mark_strategy = "samples"
    else:
        # Assume list-like of node IDs
        try:
            mark_strategy = set(int(nid) for nid in shared_nodes)
        except (TypeError, ValueError) as e:
            raise ValueError(
                "shared_nodes must be None, 'samples', or a list of node IDs"
            ) from e

    # Build TreeSequenceDictionary dict
    result = {}
    for i, (ts, idx, contig_id, symbol, typ) in enumerate(
        zip(tree_sequences, indexes, ids, symbols, types)
    ):
        # Copy and modify tables to mark shared nodes
        tables = ts.dump_tables()
        
        # Modify flags by creating a new array with selective modifications
        flags = tables.nodes.flags.copy()

        if mark_strategy == "samples":
            # Mark all sample nodes as IS_SHARED
            sample_node_ids = ts.samples()
            flags[sample_node_ids] = flags[sample_node_ids] | NODE_IS_SHARED
        elif mark_strategy is not None:
            # Mark specified node IDs as IS_SHARED (if they exist in this tree sequence)
            node_ids_to_mark = np.array(
                [nid for nid in mark_strategy if nid < ts.num_nodes],
                dtype=int
            )
            if len(node_ids_to_mark) > 0:
                flags[node_ids_to_mark] = flags[node_ids_to_mark] | NODE_IS_SHARED
        # else: mark_strategy is None, don't mark anything
        
        # Update the tables with the complete modified flags array
        tables.nodes.flags = flags

        # Create ContigKey and extract contig metadata dict
        key = ContigKey(index=idx, id=contig_id, symbol=symbol, type=typ)
        contig_dict = key._asdict()
        
        # Get current metadata
        meta = ts.metadata
        if not isinstance(meta, dict):
            meta = {}
        else:
            meta = dict(meta)  # Copy to avoid mutating original
        
        meta[CONTIG_METADATA_KEY] = contig_dict
        
        # Ensure top-level metadata schema is permissive JSON
        existing_schema = ts.metadata_schema.schema
        if not existing_schema or existing_schema.get("codec") != "json":
            tables.metadata_schema = tskit.MetadataSchema.permissive_json()
        else:
            # Ensure the schema allows the contig key
            schema_dict = copy.deepcopy(existing_schema)
            schema_dict.setdefault("additionalProperties", True)
            if "properties" not in schema_dict:
                schema_dict["properties"] = {}
            schema_dict["properties"][CONTIG_METADATA_KEY] = {
                "type": "object",
                "additionalProperties": True,
            }
            tables.metadata_schema = tskit.MetadataSchema(schema_dict)
        
        tables.metadata = meta
        
        ts_final = tables.tree_sequence()
        result[key] = ts_final

    return TreeSequenceDictionary(result)


def _get_shared_for_contig(ts, is_shared_global, contig_index):
    """
    Return a boolean array indicating which nodes are IS_SHARED and
    non-vacant for the given contig.
    
    Checks the is_vacant bitfield in node metadata to determine if
    a shared node is present (non-vacant) in this contig.
    """
    result = is_shared_global.copy()
    
    # For each shared node, check if it's vacant in this contig
    for nid in np.where(is_shared_global)[0]:
        node_meta = ts.tables.nodes[int(nid)].metadata
        if isinstance(node_meta, dict):
            vacant_bits = node_meta.get("is_vacant", 0)
            if vacant_bits & (1 << contig_index):
                result[nid] = False
    
    return result


def _try_add_vacant_flags(metadata_bytes, vacant_bitflags, schema):
    """
    Try to add ``is_vacant`` to node metadata if the schema supports it.

    Returns the (possibly updated) metadata bytes.
    """
    if vacant_bitflags == 0:
        return metadata_bytes

    schema_dict = schema.schema if schema else None
    if not schema_dict or schema_dict.get("codec") != "json":
        return metadata_bytes

    props = schema_dict.get("properties", {})
    if "is_vacant" not in props:
        return metadata_bytes

    # Decode, add field, re-encode
    try:
        if isinstance(metadata_bytes, dict):
            meta = dict(metadata_bytes)
            meta["is_vacant"] = int(vacant_bitflags)
            return meta

        if metadata_bytes:
            meta = schema.decode_row(metadata_bytes)
        else:
            meta = {}
        meta["is_vacant"] = int(vacant_bitflags)
        return schema.encode_row(meta)
    except Exception:
        return metadata_bytes

