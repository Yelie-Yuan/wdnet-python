PYTHON = python3
PACKAGE = wdnet

.PHONY: build
build:
	# pip install --upgrade build
	$(PYTHON) setup.py sdist bdist_wheel

.PHONY: clean
clean:
	rm -rf build dist src/*.egg-info src/wdnet/_*.c
	pip uninstall -y $(PACKAGE)

.PHONY: install
install:
	pip install .

venv:
	$(PYTHON) -m venv venv

# .PHONY: actvenv
# actvenv: .venv
# 	. .venv/bin/activate
