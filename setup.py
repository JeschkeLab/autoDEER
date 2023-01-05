from setuptools import setup

extras={
    "Bruker": ["XeprAPI"],
    "Matlab": ["matlabengine"],
    "Docs": ["sphinx", "pydata-sphinx-theme", "sphinx-gallery"]
}
extras["Dev"] = extras["Bruker"] + extras["Matlab"] + extras["Docs"] + \
    ["pytest"]

setup(
    name='autoDEER',
    version='0.0.1',
    author='Hugo Karas, Gunnar Jeschke and other contributors',
    package_dir={'autoDEER': 'autoDeer'},
    url = "https://github.com/HKaras/autoDeer",
    python_requires=">=3.8, <3.11",
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'deerlab >= 1.0',
        'h5py',
        'pyyaml'
    ],
    extras_require=extras
)
