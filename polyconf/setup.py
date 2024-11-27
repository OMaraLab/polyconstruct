from setuptools import find_packages, setup

setup(
    name="polyconf",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "mdanalysis",
        "click",
        "py3Dmol",
        "rdkit",
        "networkx",
        "pandas",
        "tqdm"
    ],
    python_requires=">=3.10",
)


