# Building the documentation
The autoDeer documentation is built using Spinx and the pydata theme. Please follow these instructions to build the docs.


## Installing packages

1) Install Python
2) Install autoDeer from source
3) Install Sphinx and extensions

        pip install sphinx sphinx_toolbox sphinx_copybutton numpydoc
        <!-- pip install pydata-sphinx-theme -->
        pip install furo
        pip install sphinx-gallery
        pip install sphinx-design
        pip install --upgrade myst-parser
        

## Building docs
To build docs please run:

        make clean
        make html