Installation of astroMUGS
************

How to obtain astroMUGS
=================

astroMUGS (new name of chemdiskpy) can be obtained either by cloning its Github repository or via pip.
From a terminal, type::
    
      pip install -i https://test.pypi.org/simple/ astromugs

This will install the latest version. You can use astromugs from any directory. To install the package from the Github repository, go to a directory where you want to install the code, and type:: 


    git clone https://github.com/sachagavino/astromugs.git


This will create a folder ``astromugs/``, which contains the full git repository. You can now access the package::


    cd astromugs/


To make sure you always use the latest version, you can type:: 


    git pull



Requirements and environment
=================

It is usually recommended to use a dedicated virtual environment to avoid conflicts with other packages, although astromugs requirements are widely used latest released libraries. You can use ``conda`` to create a virtual environment. 
The easiest way to do this is to use the provided ``environment.yml`` file. The name of the environment is ``astromugs`` by default, but you can change it in the ``environment.yml`` file before creating the environment.
From the terminal, type::

    conda env create -f environment.yml

Verify that the environment is created properly::

    conda env list

If the name of the environment appears in the list, it means the environment has been created successfully. 

You can now activate the new environment with::

    conda activate astromugs


Running the code
=================

In a python script, you can import 

