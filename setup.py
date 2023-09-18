from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import os
import glob

ext_data = {
    "wdnet._fib": {
        "sources": [
            os.path.join("src/_wdnet", "fib.cpp"),
            os.path.join("src/wdnet", "_fib.pyx"),
        ],
        "language": "c++"
    },
    "wdnet._utils": {
        "sources": [
            os.path.join("src/_wdnet", "utils.cpp"),
            os.path.join("src/wdnet", "_utils.pyx"),
        ],
        "language": "c++"
    },
}

extensions = []

for name, data in ext_data.items():
    sources = data["sources"]
    language = data["language"]
    obj = Extension(
        name,
        sources=sources,
        language=language,
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
    ],
    setup_requires=[
        "cython",
    ],
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
)
