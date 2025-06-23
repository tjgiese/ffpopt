#
# ------------------------------------------------------------------
# Dependency installation from conda
# ------------------------------------------------------------------
#
# 1. Download the miniforge installer
#   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
# 
# 2. Install miniforge
#   bash ./Miniforge3-Linux-x86_64.sh -b -f -p ${PWD}/miniforge3
#
# 3. Enter the conda environment
#   source ${PWD}/miniforge3/bin/activate
#
# 4. Install packages
#   conda install -y dpdata deepmd-kit ase parmed dacase::ambertools-dac=25
#
# NOTE: The precompiled dacase::ambertools-dac=25 package is not binary compatible
# with older versions of glibc that may be present on many compute clusters.
# In this case, you should compile the latest AmberTools from source.
# Alternatively, if you do not want or need to perform MM calculations with sander,
# then you can exclude AmberTools from the installation. (One still needs parmed.)
#
# For ab initio support with psi4, one can include the "psi4" package in the
# conda installation command; however, the precompiled binaries may not be
# compatible with the torch libraries installed below. In this case, you should
# either install psi4 from source, or create a separate conda environment for
# running ab initio calculations.
#
#   [See https://pytorch.org/get-started/locally/ for cuda-version-specific installation options]
#   [Currently, this appears to work for cuda 12.5 and cuda 12.6]
#
#   python3 -m pip install torch torchvision torchaudio
#   python3 -m pip install cmake tblite mace-torch geometric aimnet torchani
#
# 5. Re-enter the environment
#   conda deactivate
#   source ${PWD}/miniforge3/bin/activate
#
# 6. Use the installed software; exit with "conda deactivate"
#
#
# Note we previously used the xtb and xtb-python packages from conda; however,
# the xtb-python package is no longer maintained, and the authors recommend
# using tblite instead.
#

#
# ------------------------------------------------------------------
# Dependency installation from pip
# ------------------------------------------------------------------
#
# We will need to install 2 versions: one for pytorch and the other for tensorflow
#
# 1. Download the miniforge installer
#   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
#
# 2. Install the tensorflow python version (make sure you have gcc 9.3 (or newer) loaded to compile deepmd-kit with -std=c++20
#   bash ./Miniforge3-Linux-x86_64.sh -b -f -p ${PWD}/tensorflow
#   source ${PWD}/tensorflow/bin/activate
#   python3 -m pip install cmake==3.30 tensorflow dpdata deepmd-kit ase parmed geometric tblite
#   conda deactivate
#
# 3. Install the pytorch version (with incompatible nvidia library dependencies)
#   bash ./Miniforge3-Linux-x86_64.sh -b -f -p ${PWD}/pytorch
#   source ${PWD}/pytorch/bin/activate
#   python3 -m pip install cmake==3.30 ase parmed torch torchvision torchaudio tblite mace-torch geometric aimnet torchani
#   conda deactivate
#
# 5. Install ambertools from source
#
# 6. Install psi4
#   bash ./Miniforge3-Linux-x86_64.sh -b -f -p ${PWD}/psi4
#   source ${PWD}/psi4/bin/activate
#   conda install -y psi4
#   conda deactivate
#
# Make sure the PYTHONPATH  within ${PWD}/psi4 does not prepend those in ${PWD}/tensorflow nor ${PWD}/pytorch
# because parmed/ambertools require numpy<2, whereas conda installation of psi4 provides numpy==2
# But psi4 can still run with numpy<2
#
# If you install deepmd-kit and tensorflow via pip, then 
#


