"""
Conversion functions between a :class:`TreesArchive` and a single
:class:`tskit.TreeSequence`.

Overview
--------
**to_tree_sequence** merges all contigs into a single tree sequence by
concatenating their genomic coordinate systems, placing contigs end-to-end.
Shared nodes (IS_SHARED flag) retain their IDs; non-shared nodes get new IDs.

**from_tree_sequence** reverses this process, splitting the single tree
sequence back into a TreesArchive using the top-level metadata array that
records per-contig ``sequence_length``, ``num_nodes``, and ``contig``.

**from_slim** converts a set of SLiM-style tree sequences (where all nodes
are implicitly shared) into a :class:`TreesArchive`.
"""

import warnings

import numpy as np
import tskit

from .core import ContigKey, TreesArchive, make_contig_key, make_permissive_contig_schema
from .flags import CONTIG_METADATA_KEY, NODE_IS_SHARED

# Top-level metadata key used in the single-TS representation
_ARCHIVE_META_KEY = "tskit_multichrom_contigs"
_NULL = tskit.NULL  # -1


def to_tree_sequence(archive, record_provenance=True):
    """
    Merge a :class:`TreesArchive` into a single :class:`tskit.TreeSequence`.

    Contigs are placed end-to-end in order of their ``index`` values.
    Shared nodes (those with :data:`~tskit_multichrom.flags.NODE_IS_SHARED` set)
    retain their original node IDs. Non-shared nodes are appended after the
    last shared node for each contig and receive new IDs.

    Site and mutation metadata schemas must be identical across all contigs,
    otherwise a :class:`ValueError` is raised.

    Parameters
    ----------
    archive : TreesArchive
    record_provenance : bool
        Whether to record provenance in the output tree sequence.

    Returns
    -------
    tskit.TreeSequence
    """
    if archive.num_contigs == 0:
        raise ValueError("Cannot convert an empty archive to a tree sequence")

    sorted_keys = archive.contigs  # sorted by index

    # ------------------------------------------------------------------
    # 1. Validate that site and mutation schemas are identical
    # ------------------------------------------------------------------
    ref_site_schema = archive[sorted_keys[0]].tables.sites.metadata_schema
    ref_mut_schema = archive[sorted_keys[0]].tables.mutations.metadata_schema
    for key in sorted_keys[1:]:
        ts = archive[key]
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
    # Collect all node IDs that have IS_SHARED in at least one TS
    shared_node_ids = set()
    for key in sorted_keys:
        ts = archive[key]
        for nid in range(ts.num_nodes):
            if ts.tables.nodes[nid].flags & NODE_IS_SHARED:
                shared_node_ids.add(nid)

    # For each shared node ID, build vacant_bitflags and find first non-vacant contig
    shared_node_array = []  # list of (node_id, first_contig_key, vacant_bitflags)
    for nid in sorted(shared_node_ids):
        # vacant_bitflags: bit = contig 'index' value
        vacant_bitflags = 0
        first_key = None
        for key in sorted_keys:
            ts = archive[key]
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
    ref_ts = archive[sorted_keys[0]]
    new_tc = tskit.TableCollection(sequence_length=archive.total_sequence_length)

    # Copy individual, population tables from reference
    new_tc.individuals.replace_with(ref_ts.tables.individuals)
    new_tc.populations.replace_with(ref_ts.tables.populations)

    # Copy provenance from reference (first in index order)
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
        ts = archive[key]
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

    # ------------------------------------------------------------------
    # 5. Iterate through contigs and build node/edge/site/mutation tables
    # ------------------------------------------------------------------
    start_position = 0.0
    shared_node_counter = 0

    # node_map[contig_index][old_node_id] = new_node_id
    all_node_maps = []

    for contig_idx, key in enumerate(sorted_keys):
        ts = archive[key]
        node_map = np.full(ts.num_nodes, _NULL, dtype=np.int64)

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
                src_ts = archive[first_key]
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
                mutations.parent == _NULL,
                _NULL,
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
                tskit.NodeTableRow(flags=0, time=0.0, population=_NULL, individual=_NULL)
            )
        src_ts = archive[first_key]
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


def from_tree_sequence(ts):
    """
    Convert a single merged :class:`tskit.TreeSequence` back into a
    :class:`TreesArchive`.

    The input tree sequence must have been created by :func:`to_tree_sequence`
    (or equivalent), with top-level metadata containing a
    ``'tskit_multichrom_contigs'`` array.

    Parameters
    ----------
    ts : tskit.TreeSequence

    Returns
    -------
    TreesArchive
    """
    meta = ts.metadata
    if not isinstance(meta, dict) or _ARCHIVE_META_KEY not in meta:
        raise ValueError(
            f"Tree sequence does not have '{_ARCHIVE_META_KEY}' in top-level metadata. "
            "Was it created by tskit_multichrom.to_tree_sequence()?"
        )

    contigs_meta = meta[_ARCHIVE_META_KEY]

    # Validate cumulative lengths
    sequence_lengths = [c["sequence_length"] for c in contigs_meta]
    num_nodes_list = [c["num_nodes"] for c in contigs_meta]
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
    total_nodes = ts.num_nodes
    node_used = np.zeros(total_nodes, dtype=bool)

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

        # Determine which nodes to keep for this contig
        # Start with all nodes that have been used so far + IS_SHARED (not vacant)
        is_shared_for_contig = _get_shared_for_contig(
            ts, is_shared_global, contig_meta_dict.get("index", i)
        )
        candidate = node_used | is_shared_for_contig
        # Take the first n_nodes True values as node IDs
        candidate_ids = np.where(candidate)[0]
        if len(candidate_ids) < n_nodes:
            # Need more nodes — take next available
            all_ids = np.arange(total_nodes)
            not_used = ~candidate
            extra = all_ids[not_used][: n_nodes - len(candidate_ids)]
            candidate_ids = np.sort(np.concatenate([candidate_ids, extra]))
        else:
            candidate_ids = candidate_ids[:n_nodes]

        node_used[candidate_ids] = True

        # Extract contig tree sequence using keep_intervals then shift
        contig_ts = (
            ts.keep_intervals([[left, right]], simplify=False)
            .shift(-left, sequence_length=seq_len, record_provenance=False)
        )

        # Subset to the correct nodes
        # We need to keep candidate_ids in the new (shifted) tree sequence
        contig_tc = contig_ts.dump_tables()
        contig_tc.subset(candidate_ids)

        # Set per-contig metadata schema and metadata
        contig_schema = make_permissive_contig_schema()
        contig_tc.metadata_schema = contig_schema
        contig_tc.metadata = {CONTIG_METADATA_KEY: contig_meta_dict}

        result[key] = contig_tc.tree_sequence()

    return TreesArchive(result)


def from_slim(tree_sequences, *, slim_metadata_key="SLiM"):
    """
    Convert SLiM-style tree sequences into a :class:`TreesArchive`.

    In SLiM tree sequences all nodes are shared across chromosomes, so
    :data:`~tskit_multichrom.flags.NODE_IS_SHARED` is set on every node.
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
    TreesArchive
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
            tables.metadata_schema = make_permissive_contig_schema()
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

    return TreesArchive(result)


def _get_shared_for_contig(ts, is_shared_global, contig_index):
    """
    Return a boolean array indicating which nodes are IS_SHARED and
    non-vacant for the given contig.
    """
    # Try to read is_vacant from node metadata
    # is_vacant is a bitfield where bit i means vacant in contig i
    result = is_shared_global.copy()
    node_schema = ts.tables.nodes.metadata_schema.schema

    has_vacant = (
        node_schema
        and node_schema.get("codec") == "json"
        and "is_vacant" in node_schema.get("properties", {})
    )

    if has_vacant:
        for nid in np.where(is_shared_global)[0]:
            node_meta = ts.tables.nodes[int(nid)].metadata
            if isinstance(node_meta, dict):
                vacant_bits = node_meta.get("is_vacant", 0)
                if vacant_bits & (1 << contig_index):
                    result[nid] = False
    else:
        if np.any(is_shared_global):
            warnings.warn(
                "Node metadata does not have 'is_vacant' field; "
                "treating all IS_SHARED nodes as non-vacant for all contigs.",
                stacklevel=3,
            )

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
        if metadata_bytes:
            meta = schema.decode_row(metadata_bytes)
        else:
            meta = {}
        meta["is_vacant"] = int(vacant_bitflags)
        return schema.encode_row(meta)
    except Exception:
        return metadata_bytes
