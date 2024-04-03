from setuptools import setup, find_packages

exec(open('autodeer/_version.py').read())

extras={
    "Matlab": ["matlabengine"],
    "Docs": ["sphinx", "furo", "sphinx-gallery", "sphinx-design","myst-parser","sphinx-copybutton","sphinx-toolbox","sphinx-autoapi","sphinxcontrib-bibtex","numpydoc"],
    "GUI": ["PyQt6","threadpoolctl", "pyinstaller"],
    "test": ["pytest", "pytest-cov", "pytest-qt", "pytest-xdist"],
}
extras["Dev"] = extras["Docs"] + extras["test"]
setup(
    name='autoDEER',
    version=str(__version__),
    author='Hugo Karas, Gunnar Jeschke, Stefan Stoll and other contributors',
    package_dir={'autodeer': 'autodeer'},
    # packages=['autodeer','autodeer'],
    packages=find_packages(),
    url = "https://github.com/HKaras/autoDeer",
    python_requires=">=3.8",
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'deerlab >= 1.0',
        'pyyaml',
        'reportlab',
        'svglib',
        'xarray',
        'h5netcdf'
    ],
    extras_require=extras
)
