from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import os

ext_data = {
    "wdnet.fib_cy": {
        "sources": [
            os.path.join("src/_wdnet", "fib.cpp"),
            os.path.join("src/wdnet", "fib_cy.pyx"),
        ],
        "language": "c++"
    },
    "wdnet.utils_cy": {
        "sources": [
            os.path.join("src/_wdnet", "utils.cpp"),
            os.path.join("src/wdnet", "utils_cy.pyx"),
        ],
        "language": "c++"
    },
}

extensions = []

for name, data in ext_data.items():
    obj = Extension(
        name=name,
        sources=data["sources"],
        language=data["language"],
    )
    extensions.append(obj)

setup(
    name="wdnet",
    version="0.0.0.9000",
    author="Yelie Yuan",
    author_email="yelie.yuan@uconn.edu",
    description="A package for weighted directed networks",
    # long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Yelie-Yuan/wdnet-python",
    project_urls={
        "Bug Tracker": "https://github.com/Yelie-Yuan/wdnet-python/issues",
    },
    license="GNU General Public License (GPL)",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    # include_package_data=True,
    python_requires=">=3.10",
    install_requires=[
        "numpy",
        "pandas",
        "igraph",
        "cvxpy",
    ],
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
)
