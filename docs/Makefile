SPHINXOPTS    ?=
BUILDDIR      ?= _build

.PHONY: help html clean

help:
	@echo "Usage: make <target>"
	@echo "  html   Build HTML documentation with JupyterBook"
	@echo "  clean  Remove the build directory"

html:
	jupyter-book build . $(SPHINXOPTS)

clean:
	jupyter-book clean .
