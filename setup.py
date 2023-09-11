from setuptools import setup

exec(open('autodeer/_version.py').read())

extras={
    "Bruker": ["XeprAPI"],
    "Matlab": ["matlabengine"],
    "Docs": ["sphinx", "pydata-sphinx-theme", "sphinx-gallery", "sphinx-design"]
}
extras["Dev"] = extras["Bruker"] + extras["Matlab"] + extras["Docs"] + \
    ["pytest"]

setup(
    name='autoDEER',
    version=str(__version__),
    author='Hugo Karas, Gunnar Jeschke and other contributors',
    package_dir={'autodeer': 'autodeer'},
    packages=['autodeer','autodeer'],
    url = "https://github.com/HKaras/autoDeer",
    python_requires=">=3.8",
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'deerlab >= 1.0',
        'pyyaml',
    ],
    extras_require=extras
)
