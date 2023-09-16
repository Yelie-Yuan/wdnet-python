PYTHON = python3
PACKAGE = wdnet

.PHONY: build
build:
	# pip install --upgrade build
	$(PYTHON) setup.py sdist bdist_wheel

.PHONY: clean
clean:
	rm -rf build dist src/*.egg-info src/wdnet/_*.c src/wdnet/__pycache__
	pip uninstall -y $(PACKAGE)

.PHONY: install
install:
	pip install .
	python -c "import wdnet; \
		wdnet.hello(); \
		print(dir(wdnet)); \
		print(wdnet.fib(10));"

venv:
	$(PYTHON) -m venv venv

# .PHONY: actvenv
# actvenv: .venv
# 	. .venv/bin/activate
