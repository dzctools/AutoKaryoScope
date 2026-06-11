from setuptools import setup, find_packages

setup(
    name="AutoKaryoScope",
    version="1.0.0",
    description="AutoKaryoScope interactive chromosome-scale synteny visualization",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[],
    entry_points={
        "console_scripts": [
            "AutoKaryoScope=synteny_plot.cli:main",
        ],
    },
)
