from setuptools import setup, find_packages

setup(
    name="spare",
    version="0.1.0",
    packages=find_packages(),  # Will pick up src and rewriter as packages
    install_requires=[],
    entry_points={
        "console_scripts": [
            "compile-rewrite=rewriter.compile_rewrite_graphs:main",  # Optional CLI
        ],
    },
)