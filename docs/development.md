# Development

## Setting up a development environment

```bash
git clone https://github.com/hyanwong/tsdict.git
cd tsdict
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Running tests

```bash
pytest
```

Coverage report:

```bash
pytest --cov=tskit_multichrom --cov-report=html
```

## Code style

The project uses [ruff](https://docs.astral.sh/ruff/) for linting and
formatting:

```bash
ruff check .
ruff format .
```

## Building the documentation

```bash
pip install -e ".[docs]"
cd docs
make
```

## Contributing

Please open issues and pull requests on the
[GitHub repository](https://github.com/hyanwong/tsdict).
