PYTHON = python3
PACKAGE = wdnet

.PHONY: build
build:
	# pip install --upgrade build
	$(PYTHON) -m build

.PHONY: clean
clean:
	rm -rf build dist src/*.egg-info src/wdnet/*.c src/wdnet/__pycache__ src/wdnet/*.cpp src/wdnet/*gnu.so
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
