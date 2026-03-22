# Proposal for a tsdict library providing multiple-chromosome (contig) support for tskit

**NOTE (March 2026):** This proposal has been substantially implemented. See the "IMPLEMENTATION STATUS SUMMARY" section at the end for a detailed breakdown of what is complete and what remains to be done.

## Background

See https://docs.google.com/document/d/1mjZEuxetIfja8WD4DIt5NQXxKNe2u-kw64vjPjiwLPw

We refer to a folder of associated tree sequences as a "trees archive" (already used by SLiM), and use the generic variable name `tsd`. The object loaded from a trees archive is called a `TreeSequenceDictionary` (somewhat clunky: the name can be changed if need be). Officially each tree sequence in the archive is a "contig", but in practice, these are usually chromosomes.

Key points:

* The library is a thin Python wrapper around a dictionary of tree sequences, with the component tree sequences containing some shared data (e.g. they all have identical individuals tables)
* The dictionary values are tree sequences, with keys being namedtuples of (`index`, `id`, `symbol`, `type`) - the first 3 must be unique: `index` and `id` are integers, `symbol` and `type` are strings. However, component tree sequences can more easily be accessed via `id` or `symbol` using e.g. `ta.contig("X")`
* We aim to maintain SLiM compatibility as far as possible
* We want to support a mix of presence and absence of cross chromosomal phasing (e.g. if there are parent/child trios in a dataset, the children can be phased across chromosomes, but that may not be possible for the parents)
* Stats with `sample_sets` parameter can work across all contigs, but only if the provided sample nodes are globally phased (i.e., correspond between chromosomes); in cases where all sample nodes are globally phased in all the chromosomes, the `sample_sets` parameter can be omitted and defaults to all global samples

## Details

### Strict specifications

1. We require individual tables to be identical within an archive
2. We require population tables to be identical within an archive
3. We require migration tables to be all empty (there are better ways to record migrations now)
4. We require time_units to be identical within an archive
5. Node tables are not required to be identical (although they can be, as in SLiM), but are required to have the same metadata schema (so shared nodes can be identical, see 7 below)
6. Provenance tables are not required to be identical within an archive
7. Some nodes are shared between multiple tree sequences. We set flag IS_SHARED on nodes to specify whether this node ID is shared across all the tree sequences in this assemblage (i.e. node tables can be partially shared)
    * If one node has this flag set, each node in the archive which has both the same ID and the IS_SHARED flag set must be identical (i.e. have identical flags, times, individuals, populations, and metadata)
    * Nodes without the IS_SHARED flag are treated as specific to that local tree sequence (this includes nodes with the same ID as another IS_SHARED node on a different contig). This provides an alternative approach to the `is_vacant` bitfield flag in the SLiM representation.
8. The top level metadata of each tree sequence must contain a 'contig' key in the same format as used by the `this_chromosome` field in SLiM, that specifies "index", "id", "symbol", and "type". "Id" and "index" are required to be non-negative integers, and unique across contigs. The order of contigs in the archive is determined by the order of the "index" values, but there is a looser specification for "index" in that the integers are not required to be consecutive (see "reindexing" below).

### Optional table specifications

The following are optional but highly recommended (reasons explained below)

9. Site and mutation metadata schemas should be identical within an archive. If not, site (or mutation) tables cannot be reliably merged, so an archive cannot be converted into a single tree sequence representation.
10. (Currently) Top level metadata should be JSON format with a permissive schema. If not, then metadata will not be fully round-tripped (in particular, the original schema of each contig is not saved anywhere when converting to a single tree sequence)
11. The node metadata schema should be compatible with storing the is_vacant metadata information used by SLiM. If not, data can be converted to a single tree sequence format, but it is not guaranteed to retain the same node IDs and sample flags when being converted back to archive format.

### Constraints

The tskit node ID is used to link all nodes that are shared. This has the major advantage that a single ID is valid both for each constituent tree sequence and for the total assemblage. The down side to this is that, as node IDs correspond to their order in the node tables, having a shared node ID of (say) 99 requires tree sequences sharing that node to have at least 100 nodes in total. We can minimise this problem by using low IDs for "shared nodes" (likely to be generally true, as shared nodes are most often samples, and samples often have IDs 0..N). If necessary we can also pad the node table with unused nodes. 

## Basic functions

### Loading and caching

✅ **IMPLEMENTED**

Loading: `tsdict.load(filename)`

1. ✅ We check that the tree sequences in the archive conform to the requirements above, in particular checking for node identity requirements (point 7 above)
2. ✅ We perform caching (see below)

Caching (also performed when creating a new archive via .subset)

1. ✅ We cache the total_sequence_length of all contigs. Note that this might overcount genome size, e.g. it adds *both* X and Y lengths to the total.
2. ✅ We calculate a cache of the IDs of nodes which are cross-phased over all the contigs. These IDs are valid for use in whole-genome statistical calculations. 
3. ✅ We also count the number of "nonglobal sample nodes": i.e. nodes that are samples but not shared (cross-phased) in one or more contigs. If there are any nonglobal sample nodes, we can't easily refer to a single sample set. An archive with any nonglobal sample nodes is called a nonglobal-sample ARG. (This term is specifically about sample nodes; nonsample nodes may be non-cross-phased in any case.)
4. ✅ We cache (and index) contig metadata information so that it is easy to refer to the tree sequences corresponding to individual contigs, e.g. via `ta.contig("X")`

### Subsetting

✅ **IMPLEMENTED**

We should be able to create a new archive by subsetting. This allows us to create an object without any nonglobal sample nodes (for example, all the autosomes using `ta.subset(type="A")`). The subsetting can be done entirely in-memory, and should not require making a copy of the underlying tree sequences, so should be cheap. It should only require re-running the caching process (above).

### Reindexing

✅ **IMPLEMENTED**

Even if a set of tree sequences are indexed from 0..N, when subsetted those indexes may not remain consecutive integers. A slightly more expensive operation involves making a copy of the tree sequences and changing the indexes in their metadata to be from 0..N. We could also use a `reindex` function to reorder the contigs, by specifying a different index order.

### Conversion

✅ **MOSTLY IMPLEMENTED** — `to_ts()`, `from_ts()`, `from_slim()`, and `from_tree_sequences()` functions exist.

From a SLiM tree sequence archive:

All we really need to do is to:
1. Mark the IS_SHARED flag on all nodes
2. Duplicate or alias the top level "['SLiM']['this_chromosome']" metadata (with keys "index", "id", "symbol", and "type") into a ['contig'] metadata field stored at the top level
3. (possibly) translate the is_vacant bitfield into the partly-shared node format. This simply means that for any node which is marked vacant in a specific contig, the node in the table for that contig has:
    * The IS_SHARED flag unset
    * The individual value set to NULL
    * (probably) the sample flag removed - this means that it won't be treated as a sample for stats purposes etc.


To a single tree sequence

We only allow merging if site and mutation metadata schemas are all identical, otherwise we raise an error. Contigs are taken in the order of their `id` in the top level metadata
* Top level metadata: The metadata for each individual contig is placed as an item in a top level JSON array in the new file, with a permissive_json schema.
* Additionally the `sequence_length` and `num_nodes` for each of the original contigs is added to their metadata.
* The migrations tables is kept blank
* We make a "shared_node_array" containing, for each node that has an IS_SHARED flag in at least one of the tree sequences, a tuple of (node ID, first_chrom_id, vacant_bitflags), where each bit in vacant_bitflags corresponds to the metadata 'index' of a contig, as described in the SLiM manual. Here, first_chrom_id is the first contig index where the node is present and not vacant; absent indexes are set to 1 in vacant_bitflags. The array is sorted by node ID, and we set a variable shared_node_counter=0.
* We create a new TableCollection and fill:
   * The individual and population tables from any of the tree sequences (all are identical)
   * The provenance table from the first tree sequence (the one with the lowest "index"). This means for simplicity the new TS will only have the provenance of the first contig, and we don't guarantee to be able to round-trip the provenances.
   * A variable start_position = 0 
* We iterate through all the contigs in order of `index`.
   * Using np.arange, we make arrays mapping of all the nodes in the current contig to their current ID.
   * We iterate through the contig's node table, iteratively adding nodes from there to the target node table in the new TableCollection.
      a. If the number of rows in the target node table is equal to ID = shared_node_array[shared_node_counter][0], we find the first non-vacant  shared_node_chromosome[shared_node_counter], append the node with that ID from there, and increment the shared_node_counter. If the metadata allows, we also save the vacant_bitflags into that node metadata, under the "is_vacant" key.
      b. Otherwise we add the node from the current contig, and save the new ID as the value in the map.
   * We add the chromosome edges using `append_columns` with left=left+start_position, right=right + start position, and parent = node_map[old_parent], child=node_map[old_child].
   * We add the chromosome's mutations using `append_columns` with node=node_map[old_node], parent=np.where(parent==NULL, NULL, parent+new_tc.mutations.num_rows) and site=site+new_tc.sites.num_rows.
   * We add the chromosome's sites using `append_columns` setting position = old_position + start_position
   * We add the chromosome's sequence_length to start_position
	(note that there may be an equivalent way to do this using ts.concatenate)

Note that there may be many more sample nodes in the resulting single tree sequence than there were in any one of the original tree sequences, as non-shared sample nodes will all have a new ID (we could assert the number of these, as a check).

From a single tree sequence

We assume that the tree sequence has a top-level metadata consisting of an array, each element of which describes a chromosome, and which contains a numerical "sequence_length", "num_nodes", and a "chromosome" dictionary.

We check that the `index` in each chromosome key is increasing in order in the top level metadata array, store cumsum=[0] + [cumulative sum of the sequence_lengths of the chromosomes], and check that the last entry adds up to the sequence_length of the current single tree sequence.

We make a node_used boolean array of the length of the total number of nodes in the original single tree sequence.

We enumerate through the top-level array getting `i` and the `num_nodes` for each chromosome:
* For the edge table, perform keep_intervals(cumsum[i], cumsum[i+1]) on the original, then .shift(-cumsum[i], sequence_length-cumsum[i+1] - cumsum[i]), making sure these don't simplify(). This should return a chromosome-level tree sequence with the correct site and mutation counts.
* The individual, population, and provenance tables should be automatically copied
* The tricky thing is to know which nodes to keep. We need to make sure that all the nodes marked as IS_SHARED and which do not have an is_vacant bitflag for this chromosome in their metadata retain the same IDs in the final tree sequence. This can be done using .subset. Specifically, we make a boolean is_shared array of the same length as node_used, but full of zeros, and set True for those nodes marked as IS_SHARED and which do not have the is_vacant bitfield set (we warn if there is no is_vacant bitfield, and treat this as not being vacant). We logical_or the is_shared and node_used arrays, and take the indexes of the first num_nodes True values as the node IDs to pass to subset().

### Simplifying

✅ **PARTIALLY IMPLEMENTED**

`TreeSequenceDictionary.simplify(samples=None, *, individuals=None)` simplifies all contigs in a coordinated way.

* `samples=[...]`: provided node IDs must be globally phased.
* `individuals=[...]`: simplifies each contig using sample nodes for those individuals; this works on nonglobal-sample ARGs where nonglobal sample nodes are present.
* With neither argument, simplification defaults to global phased samples and is allowed only when the assemblage is not a nonglobal-sample ARG.
* ❌ other arguments to simplify not yet implemented (details to work out: some
arguments may not be appropriate)

## General API

✅ **PARTIALLY IMPLEMENTED**

* ✅ Stats with `sample_sets` parameter (e.g. `ta.stats.diversity(sample_sets=[[3,4,5,6]])`) can work across all contigs, but only if sample node IDs are globally phased. Basic stats like `diversity()` are implemented
  - ❌ Windowing support is not yet implemented for cross-contig stats (raises NotImplementedError)
  - ❌ Other cross-contig stats methods not yet implemented
* ✅ `ta.shared_node_ids` property provides globally shared node IDs
* ✅ `global_phased_node_ids` property lists globally phased node IDs
* ✅ Contig accessor: `ta.contig(id_or_symbol)` — access by integer id, symbol, or index
* ✅ Iterator: `for key, ts in ta.items():`
* ✅ Dictionary-like interface: `ta.keys()`, `ta.values()`, `ta[key]`
* ❌ `.variants()` iterator — not yet implemented
* ❌ `.trees()` iterator — not yet implemented
* ❌ Recapitation support (combining multi-chromosome SLiM tree sequences from parallel msprime simulations) — not yet implemented; would use `union` or similar

---

## IMPLEMENTATION STATUS SUMMARY

### Core Functionality (✅ Complete)
- ✅ TreeSequenceDictionary class with validation
- ✅ ContigKey namedtuple for organizing contigs
- ✅ Loading and saving archives (load/dump)
- ✅ Caching system (total_sequence_length, global_phased_node_ids, nonglobal_sample_node_count, contig metadata)
- ✅ Node IS_SHARED flag support
- ✅ Strict specifications enforcement (individual/population table identity, migration table emptiness, time_units, node metadata schema consistency)
- ✅ Optional table specifications (site/mutation metadata schema validation)
- ✅ Dictionary-like interface (keys(), values(), items(), __getitem__)
- ✅ Contig access by id, symbol, or index

### Subsetting and Reindexing (✅ Complete)
- ✅ Subset by symbols, type(s), ids, or indexes
- ✅ Reindex to renumber from 0..N
- ✅ Reorder contigs by specifying desired order

### Conversion Functions (✅ Complete)
- ✅ ta.to_ts() — merge TreeSequenceDictionary into single TreeSequence
- ✅ from_ts() — split single TreeSequence into TreeSequenceDictionary
- ✅ from_slim() — convert SLiM-style tree sequences to TreeSequenceDictionary
- ✅ from_tree_sequences() — create assemblage from list of tree sequences

### Statistics (⚠️ Partial)
- ✅ Cross-contig stats with `sample_sets` parameter (requires globally phased nodes) — `diversity()` implemented
- ❌ Windowing support for cross-contig stats (raises NotImplementedError)
- ❌ Other cross-contig stats methods not yet implemented

### Simplification (✅ Complete)
- ✅ ta.simplify(samples=None, *, individuals=None)
- ✅ Supports simplification of nonglobal-sample ARGs via `individuals=[...]`

### Missing/Incomplete Features (❌)
- ❌ .variants() iterator
- ❌ .trees() iterator
- ❌ Windowed statistics
- ❌ Recapitation/union support for multi-chromosome SLiM sequences

### Notes
- All strict specifications (1-8) are enforced
- Optional table specifications (9-11) are validated
- The implementation has solid test coverage in the tests/ directory
- The API matches the proposal well for core functionality
