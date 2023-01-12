from setuptools import setup

exec(open('autoDeer/_version.py').read())

extras={
    "Bruker": ["XeprAPI"],
    "Matlab": ["matlabengine"],
    "Docs": ["sphinx", "pydata-sphinx-theme", "sphinx-gallery"]
}
extras["Dev"] = extras["Bruker"] + extras["Matlab"] + extras["Docs"] + \
    ["pytest"]

setup(
    name='autoDEER',
    version=str(__version__),
    author='Hugo Karas, Gunnar Jeschke and other contributors',
    package_dir={'autodeer': 'autoDeer'},
    packages=['autodeer','autodeer'],
    url = "https://github.com/HKaras/autoDeer",
    python_requires=">=3.8, <3.11",
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'deerlab >= 1.0',
        'pyyaml'
    ],
    extras_require=extras
)
