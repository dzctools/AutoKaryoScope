from setuptools import setup, find_packages

setup(
    name="synteny_plot",
    version="1.0.0",
    description="Interactive synteny plot with SV detection",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[],
    entry_points={
        "console_scripts": [
            "synteny-plot=synteny_plot.cli:main",
        ],
    },
)
