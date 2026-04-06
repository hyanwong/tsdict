SPHINXOPTS    ?=
BUILDDIR      ?= _build
PROJECT_ROOT  := $(abspath ..)

.PHONY: help html clean

help:
	@echo "Usage: make <target>"
	@echo "  html   Build HTML documentation with JupyterBook"
	@echo "  clean  Remove the build directory"

html:
	PYTHONPATH="$(PROJECT_ROOT):$$PYTHONPATH" jupyter-book build . $(SPHINXOPTS)

clean:
	jupyter-book clean .
