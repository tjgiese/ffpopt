#
# ------------------------------------------------------------------
# Installation from conda
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
#   [See https://pytorch.org/get-started/locally/ for cuda-version-specific installation options]
#   [Currently, this appears to work for cuda 12.5 and cuda 12.6]
#
#   python3 -m pip install torch torchvision torchaudio
#   python3 -m pip install tblite mace-torch geometric aimnet torchani
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
