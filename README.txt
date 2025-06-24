# ------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------
#
# 1. Install the dependencies. Most dependencies are optional.
#    Optional packages (technically)
#      a. sander/pysander (provides --model=sander)
#      b. dpdata & deepmd-kit & tensorflow (provides --model=qdpi2)
#      c. torchani, aimnet, mace are all pytorch models providing
#         --model=mace --model=aimnet2, --model=aimnet2_wb97m 
#                      (aimnet2 is an alias for aimnet2_wb97m)
#         --model=aimnet2_b973c --model=aimnet2_qr
#         --model=ani1x --model=ani2x --model=ani1ccx
#      d. If the pytorch models are not used, then torch, 
#         torchvision, and torchaudio are optional as well.
#      e. tblite (provides --model=xtb and needed for --model=qdpi2)
#      f. psi4 (profides --model="theory/basis")
#
#  2. Install fppopt
#      cd build
#      bash ./run_cmake.sh
#
#  If you have a modern OS, it's easiest to install the dependencies
#  using conda; however, note that the psi4 conda package will install
#  numpy-2.x, whereas parmed requires numpy-1.x.  You will need to
#  create a separate conda installation specifically for psi4. If
#  you install psi4 in miniforge3_psi4 and all other dependencies in
#  miniforge3_allother, then you can *prepend* the bin and site-packages
#  from miniforge3_allother to your PATH and PYTHONPATH variables,
#  and then *append* the bin and site-packages from miniforge3_psi4.
#  The psi4 pacakge will still run with numpy-1.x even though it
#  installs numpy-2.x; and the procedure that I've described will
#  cause the numpy-1.x installation to be loaded rather than 2.x.
#
#  If you are on a cluster that has an "old" glibc version, then
#  the  dacase::ambertools-dac=25 package may not work correctly.
#  One cannot use the unofficial "ambertools" conda package because
#  that version of the code has a bug in the charmm module that
#  produces memory leaks and segmentation faults when the calculator
#  is reloaded more than 1 time.  For older machines, you should
#  install ambertools from source and install the other dependencies
#  using pip. You will need to create separate installations for
#  tensorflow and pytorch models. The tensorflow and pytorch packages
#  require different versions of the nvidia libraries that are
#  incompatible with each other.
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
#
# ------------------------------------------------------------------
# Documentation
# ------------------------------------------------------------------
#
# All scripts support a collection of common options, and each script
# includes a few script-specific options.
#
# The common options are:
#   --model : str
#       Select the model chemistry. The available options depend on
#       which dependencies were installed. The default value is
#       "sander". Other options include xtb, qdpi2, mace, aimnet2,
#       aimnet2_wb97m (aimnet2 is an alias for aimnet2_wb97m),
#       aimnet2_b973c, aimnet2_qr, ani1x, ani2x, ani1ccx, theory/basis
#
#  --no-opt
#       Calculate a single point energy rather than a geometry
#       optimization
#
#  --ase-opt
#       Use the ASE optimzer rather than geomeTRIC
#
#  --geometric-maxiter: int
#       Number of geometry optimization steps when using geometric.
#       Default: 500
#
#  --geometric-coordsys: str
#       Coordinate system. Default: tric
#
#  --geometric-converge: str
#       Optimization tolerance. Default: 'set GAU_TIGHT'
#
#  --geometric-enforce: float
#       Constraint enforcement tolerance. Default 0.1
#
#  --psi4-memory: str
#       Amount of RAM used by psi4. Default: '1gb'
#
#  --psi4-num-threads: int
#       The number of psi4 threads. Default: 4
#
#  --parm: str
#       The amber parm7 file. This is always required.
#
#  --crd: str
#       The amber formatted rst7 file. This is always required.
#
# 
# -----------------
# ffpopt-Optimze.py
# -----------------
#   Reads parm7 & rst7, performs a geometry optimization, and writes
#   the structure and energy to --oscan=oscan.xyz
#   If --no-opt is provided, then it calculates a single point energy.
#   If --iscan=iscan.xyz is provided, then it performs a geometry
#   optimization of each structure within iscan.xyz.
#   If --iscan and --no-opt are provided, then it performs a single
#   point calculation for each structure.
#   If --ase-opt is given, then the ASE BFGS optimizer is used
#   (the default is to optimize with geomeTRIC).
#
#   --oscan: str
#       The output xyz filename
#
#   --iscan: str
#       Optional. The input xyz filename. If not provided, the
#       coordinates within --crd are used.
#
#   --ignore
#       If present, then ignore the constraint definitions stored
#       within the xyz file. The default behavior is to enforce
#       the constraints specified within --iscan.
#
#   --constrain: str
#       A comma-separated list of 2, 3, or 4 zero-based integers
#       (the first atom is index 0).  The list can be appended with
#       =value to specify a constraint value, otherwise the constraint
#       value is the initial value calculated from the input
#       coordinates. This option can be used multiple times to enforce
#       multiple constraints.
#       For example,
#       --constrain='0,1=2.0' will constrain the bond between the
#       first 2 atoms to be 2.0 Angstroms.
#       --constrain='0,1,2=30.' will constrain the angle between the
#       first 3 atoms to be 30. degrees.
#       A list of 4 atoms constrains a dihedral.
#
# -------------------
# ffpopt-DihedScan.py
# -------------------
#   Reads parm7 & rst7 and a list of 4 atoms defining a dihedral angle.
#   It then scans the dihedral with a series of relaxed optimizations.
#   The procedure is to: 1. optimize for a minimum. 2. Create a schedule
#   of angles that uniformly span [0,360). 3. Find the position in the
#   schedule that best matches the optimized dihedral. 4. Sequentially
#   optimize the dihedrals in the schedule in the foward direction until
#   all 360 degrees are considered. 5. Repeat the scan in the reverse
#   direction. 6. Sort the two scan and choose the geometry & energy
#   that produced the lowest energy.  The output structures are written
#   to --oscan.
#
#   --dihed: str
#        A comma-separated list of 4 zero-based integers defining a
#        dihedral angle.
#
#   --oscan: str
#        The output XYZ file
#
#   --delta: int
#        The scan spacing. Default: 10 degrees
#
#   --constrain: str
#        Additional constraints applied to each structure.
#        These can be bonds, angles, or dihedrals. These coordinates
#        are not scanned. See ffpopt-Optimize.py for an extended
#        description. This option can be used more than once.
#
#
# ----------------------
# ffpopt-GenDihedFit.py
# ----------------------
#   Reads a json input file and optimize torsion parameters.
#
#   inp: str, positional argument
#        The name of the json file.
#
#   --stride: int
#        The stride used when reading the input scans. Default: 1
#
#   --nlmaxiter: int
#        Maximum number of nonlinear optimization steps. Default: 200
#
#   --nlrhobeg: float
#        Initial parameter displacements. Default: 0.25 kcal/mol.
#
#   --nltol: float
#        Parameter optimization termination tolerance. Default: 0.01.
#
#
# -----------------------------
# ffpopt-DihedTwistWorkflow.py
# -----------------------------
#   Reads a json input file and writes a bash script that uses
#   ffpopt-DihedScan.py and ffpopt-GenDihedFit.py to iteratively
#   parametrize torsion potentials.
#   
#   --bond: str
#        Two 0-based integers separated by a comma. This option
#        can be used more than once.
#
#   --delta: int
#        The scan spacing. Default: 10 degrees
#
#   --nprim: int
#        The number of primitive torsion functions applied to each
#        dihedral. Default: 3
#
#   --maxiter: int
#        The number of training iterations. Each iteration repeats
#        the sander scans with the current set of parameters.
#        Default: 2
#
#   --bytype
#        If present, then all parameters are based on atom types,
#        and applied globally (even if the atom quartet isn't being
#        scanned).  The global parameters are also written to a
#        frcmod file.
