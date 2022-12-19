Installation
===================

Currently, the only way to install autoDEER is from source.

Installing from source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To download the latest version using gig. ::

    git clone https://github.com/JeschkeLab/DeerLab.git

The latest version can **always** be found in the *main* branch ::

    git pull origin main

AutoDEER can eiter be installed normally or in develop mode. Develop mode is
useful if you plan to develop/expand autoDEER. ::
        
    python -m pip install .

Or to install in develop mode. ::
        
    python -m pip install -e .  


Installing Python CentOS 7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Warning::
    CentOS 7 will be officially discontinued on the 30th June 2024, after this
    date there is an increased security risk with the continued use of this Operating
    System. 

.. Caution::
    Installing Python requires root privileges. 

On Bruker spectrometers the computer is closely tied to the console, for this
reason it is often the case that OS on the computer becomes out of date. One 
common, OS to be stuck with is CentOS 7. CentOS was, untill 2021, the free
version of Red Hat Enterprise Linux (RHEL), which in turn is a long term support
(lts) version of the eminent linux dstribution Fedora. CentOS 7 was released
in 2014 based on Fedora 19, released in 2013. 

The problem
    Linux OSs come with python pre installed, or to be specific Python2 installed.
    This is because python is a necessary dependecy for the linux OS, however we require
    much more modern versions of python. Installing such versions can be a challenge. 

The solution
******************

.. danger:: 
    This solution requires the installation of software outside the sphere
    of the development team. The developer of autoDEER take no responsibility
    for any issues caused by these instructions. A full system backup is recommended
    before attempting this install.

The first step to installing anything on a Linux OS, or any operating system, 
is to perfom a system and repository update. ::

    yum -y update

The most fundamental issue with installing python on CentOS is a conflict 
between the openSSL versions. openSSL is an IP encryption package.::

    yum install make gcc perl pcre-devel zlib-devel
    wget https://ftp.openssl.org/source/old/1.1.1/openssl-1.1.1.tar.gz
    tar xvf openssl-1.1.1.tar.gz
    cd openssl-1.1.1/
    ./config --prefix=/usr --openssldir=/etc/ssl --libdir=lib no-shared zlib-dynamic
    make
    make test
    make install

Once this is complete we need to add the new openSSL to our path.::
   
    export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64
    echo "export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64" >> ~/.bashrc

To check that it has installed correctly type: ::
    
    openssl version

This should now be `OpenSSL 1.1.1`

Before we can install python we will want to install some extra dependecies,
these are required before XeprAPI can be insatlled. ::
   
    yum install tcl
    yum install tcl-devel

    yum install tk
    yum install tk-devel

Now we have the necessary dependecies for installing python. However, we need 
a system to manage multiple python versions. The recommended way of doing this
is using pyenv (python enviroments).


Now that pyenv has been installed we can install the python version that we want.
In this case it is 3.9.7::

    pyenv install 3.9.7

We then want to set the newly installed pyenv to be default for our directory::

    pyenv local 3.9.7


FAQs
******************
1. Out of date SSL/CA certificates:
    sudo ./usr/sbin/update-ca-certificates
