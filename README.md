# InsFEM
Welcome to the world of data analysis and problem reduction by means of revolutionary techniques and code!

## INSTALL
Currently InspiraFEM uses the following compilers:
GCC-9.0

Currently InspiraFEM uses the following packages:

**MPICH**:
http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz

**PETSc**:
http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.13.3.tar.gz

**SLEPc**:
https://slepc.upv.es/download/distrib/slepc-3.13.3.tar.gz

**VTK**:
https://www.vtk.org/files/release/9.0/VTK-9.0.1.tar.gz

**HDF5**:
https://www.hdfgroup.org/package/hdf5-1-12-0-tar/?wpdmdl=14581&refresh=5f19d2bd301851595527869

**CONDA or MINICONDA**:
We recommend making your own conda environment for InspiraFEM,

conda create --name inspirafem

Locate the directory for the conda environment in your terminal window. 
Enter that directory and create these subdirectories and files: 

    cd $CONDA_PREFIX
    mkdir -p ./etc/conda/activate.d
    mkdir -p ./etc/conda/deactivate.d
    touch ./etc/conda/activate.d/env_vars.sh
    touch ./etc/conda/deactivate.d/env_vars.sh

You can also use the ones included in [activate](https://github.com/InspiraSM/temporal_repo/blob/master/activate.d/env_vars.sh) 
and [deactivate](https://github.com/InspiraSM/temporal_repo/blob/master/deactivate.d/env_vars.sh)
This way conda will manage your environment variables without mangling the rest of your paths

You will additionally need the following Perl libraries that can be installed with cpanm

    sudo apt install cpanminus

    sudo cpanm install File::Find::Rule
    sudo cpanm install File::Copy::Recursive
    sudo cpanm install Parallel::ForkManager
    sudo cpanm install Sort::Naturally
    sudo cpanm Recursive

You may need to be sudo for the prevous installations.


After all of this is installed git clone into the InspiraFEM repository and execute

make (which will implicitly make config, debug,release and test, you can execute each of this separately)

After the tests run (and pass) you are ready to go!
