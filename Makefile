PYTHON = python3

.PHONY: build
build:
	pip install --upgrade build
	$(PYTHON) -m build

.PHONY: clean
clean:
	rm -rf dist *.egg-info

.PHONY: install
install:
	pip install --editable .

venv:
	$(PYTHON) -m venv venv

# .PHONY: actvenv
# actvenv: .venv
# 	. .venv/bin/activate
