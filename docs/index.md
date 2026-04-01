# tskit_multichrom

A Python library for efficiently storing and analysing multiple chromosomes
(contigs) using [tskit](https://tskit.dev).

The central object is {class}`~tskit_multichrom.TreeSequenceDictionary`
(abbreviated `tsd`), which holds a collection of
{class}`tskit.TreeSequence` objects — one per contig — in a dictionary keyed
by {class}`~tskit_multichrom.ContigKey` named-tuples.  It provides:

- **Unified storage** — manage all chromosomes as a single object.
- **Cross-chromosome statistics** — compute diversity, PCA, etc. over shared
  sample nodes.
- **Efficient subsetting** — filter by chromosome type without copying data.
- **Format conversion** — round-trip between multi-contig archives and a
  single merged tree sequence.
- **SLiM compatibility** — import SLiM per-chromosome tree sequences directly.

::::{grid} 1 1 2 2

:::{grid-item-card} Getting started
:link: installation
:link-type: doc

Install the library and follow the quick-start guide.
:::

:::{grid-item-card} API reference
:link: api
:link-type: doc

Full documentation for every public class and function.
:::

::::
