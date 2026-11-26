"""Compile the private Cython binding and portable RNAlib sources."""

from pathlib import Path

from Cython.Build import cythonize
from setuptools import Extension, setup

vendor = Path("vendor/viennarna")
sources = (vendor / "sources.txt").read_text().splitlines()
setup(
    ext_modules=cythonize(
        [Extension(
            "DeepOneCalculator._viennarna",
            ["src/DeepOneCalculator/_viennarna.pyx", *sources],
            include_dirs=[str(vendor), str(vendor / "src"), str(vendor / "src/ViennaRNA")],
            define_macros=[("HAVE_CONFIG_H", "1"), ("VRNA_LOG_NO_DEBUG_RNALIB", "1")],
            libraries=["m", "pthread"],
            extra_compile_args=["-O2", "-fvisibility=hidden"],
        )],
        compiler_directives={"language_level": 3},
    ),
)
