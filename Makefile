PYTHON = python3
PACKAGE = wdnet

.PHONY: build
build:
	# pip install --upgrade build
	$(PYTHON) setup.py sdist bdist_wheel

.PHONY: clean
clean:
	rm -rf build dist src/*.egg-info src/wdnet/_*.c src/wdnet/__pycache__ src/wdnet/*.cpp
	pip uninstall -y $(PACKAGE)

.PHONY: install
install:
	pip install .
	python -c "import wdnet; \
		wdnet.hello(); \
		print(dir(wdnet)); \
		print(wdnet.fib(10));\
		from wdnet._utils import node_strength_py as s;\
		print(s([(1, 2), (2, 3)], [0.5, 10]));\
		from wdnet import WDNet;\
		WDNet(edgelist=[(1, 2), (2, 3)], edgeweight=[2, 10])"

venv:
	$(PYTHON) -m venv venv

# .PHONY: actvenv
# actvenv: .venv
# 	. .venv/bin/activate
