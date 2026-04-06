.. _sec_api:

=================
API Documentation
=================

.. currentmodule:: tsgroup


.. _sec_api_toplevel:

*****************************
Top-level functions
*****************************

.. autofunction:: load

.. autofunction:: from_tree_sequences

.. autofunction:: from_ts

.. autofunction:: from_slim


.. _sec_api_tsd:

****************************
TreeSequenceGroup
****************************

.. autoclass:: TreeSequenceGroup
    :members: validate, contig, subset, reindex, simplify, dump, to_ts,
              keys, values, items,
              contigs, num_contigs, total_sequence_length,
              global_phased_node_ids, shared_node_ids,
              nonglobal_sample_node_count, is_nonglobal_sample_arg,
              stats


.. _sec_api_stats:

****************************
Statistics
****************************

The statistics sub-namespace is accessed via
:attr:`TreeSequenceGroup.stats`:

.. autoclass:: tsgroup.stats.TreeSequenceGroupStats
    :members: diversity, pca


.. _sec_api_contig_key:

****************************
ContigKey
****************************

.. autofunction:: make_contig_key

.. autoclass:: ContigKey


.. _sec_api_flags:

****************************
Flags and constants
****************************

.. autodata:: NODE_IS_SHARED

.. autodata:: CONTIG_METADATA_KEY

.. autodata:: ARCHIVE_EXTENSION
