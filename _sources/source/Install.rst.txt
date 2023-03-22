Installation
===================

Currently, the only way to install autoDEER is from source.

Installing from source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To download the latest version using gig. ::

    git clone https://github.com/JeschkeLab/DeerLab.git

The latest version can **always** be found in the *main* branch ::

    git pull origin main

Multiple versions of autoDEER exist depending on what hardware you use:
    
1. Bruker
2. Matlab
3. Docs
4. Dev (Install all of the above plus pytest)

To install the Bruker module, for example::

    python -m pip install .[Bruker]

When using MacOS or any ZSH based shell you need to do::

        python -m pip install ".[Bruker]"

AutoDEER can also be installed in editable mode. Editable mode is
useful if you plan to develop/expand autoDEER. ::

    python -m pip install -e .  


.. toctree::
    .. toctree::
    :hidden:
    :maxdepth: 2
    
    ./Install_python.rst