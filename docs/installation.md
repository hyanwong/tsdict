# Installation

## From PyPI

```bash
pip install tskit_multichrom
```

## From source

```bash
git clone https://github.com/hyanwong/tsdict.git
cd tsdict
pip install -e ".[dev]"
```

## Building the documentation

Install the docs extras and then run JupyterBook from the `docs/` directory:

```bash
pip install -e ".[docs]"
cd docs
make
```

The built HTML will appear in `docs/_build/html/`.
