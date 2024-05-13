from setuptools import find_packages, setup

setup(
    name="polytop",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "mdanalysis",
        "click",
        "py3Dmol",
        "rdkit",
    ],
    python_requires=">=3.10",
)

if __name__ == "__main__":
    setup()
