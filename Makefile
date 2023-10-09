PYTHON = python3
PACKAGE = wdnet

.PHONY: build
build:
	# pip install --upgrade build
	$(PYTHON) -m build

.PHONY: clean
clean:
	rm -rf build dist src/*.egg-info src/wdnet/*.cpp src/wdnet/*/*.cpp src/wdnet/__pycache__
	pip uninstall -y $(PACKAGE)

.PHONY: install
install:
	pip install .
	$(PYTHON) test_on_install.py

venv:
	$(PYTHON) -m venv venv

# .PHONY: actvenv
# actvenv: .venv
# 	. .venv/bin/activate
