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
#  2. Install ffpopt
#      cd build
#      bash ./run_cmake.sh
#      (by default, run_cmake.sh installs to ../local)
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
#   conda install -y dpdata deepmd-kit ase parmed dacase::ambertools-dac=25 dftd3-python
#
# NOTE: The precompiled dacase::ambertools-dac=25 package is not
# binary compatible with older versions of glibc that may be present
# on many compute clusters. In this case, you should compile the
# latest AmberTools from source. Alternatively, if you do not want
# or need to perform MM calculations with sander, then you can exclude
# AmberTools from the installation. (One still needs parmed.)
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
# We will need to install 2 versions: one for pytorch and the
# other for tensorflow
#
# 1. Download the miniforge installer
#   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
#
# 2. Install the tensorflow python version (make sure you have
#    gcc 9.3 (or newer) loaded to compile deepmd-kit with -std=c++20
#   bash ./Miniforge3-Linux-x86_64.sh -b -f -p ${PWD}/tensorflow
#   source ${PWD}/tensorflow/bin/activate
#   python3 -m pip install cmake==3.30 tensorflow dpdata deepmd-kit ase parmed geometric tblite dftd3
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
#   conda install -y psi4 simple-dftd3 dftd3-python
#   conda deactivate
#
# Make sure the PYTHONPATH  within ${PWD}/psi4 does not
# prepend those in ${PWD}/tensorflow nor ${PWD}/pytorch
# because parmed/ambertools require numpy<2, whereas conda
# installation of psi4 provides numpy==2 But psi4 can still
# run with numpy<2
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
#
#
# -------------------------------------------------
# JSON input file format for ffpopt-GenDihedFit.py
# -------------------------------------------------
# The json file contains 3 main keys: params, output, and systems.
#
# The params dictionary provides a petite list of unique dihedral
# parameters.  Its subkeys are "nprim" and "masks".
#
# The nprim value is the integer number of primitive dihedrals.
# For example, nprim=3 would model the dihedral with 3 torsion
# potentials with periodicities 1, 2, and 3. The corresponding
# 3 force constants are to be determined.
#
# The masks value is either "null" or a list containing
# sublists of 4 atom type amber masks. If the value is null,
# then the parameter is
# "bespoke"; that is, the potential must be manually mapped
# to atom name quartets. In constrast, if a list of atom type
# masks are provided, then the potential is applied to all
# proper torsions with a matching set of atom types -- even
# if those torsions are not being scanned.  These are "global"
# parameters that can be written into a frcmod file and applied
# via tleap.
# An example, to define a global parameter applied to all
# cd-nf-ce-o torsions, one would set "masks": [ ["@%cd","@%nf","@%ce","@%o"] ].
# The value of masks is a list of lists because one could choose
# to parametrize a single potential for multiple atom type
# quartets.
#
# The "output" value is the name of an output frcmod file.
# All global parameters will be written to the frcmod file.
# If all parameters are bespoke, then the frcmod file is empty.
#
# The "systems" value is a list of dictionaries. Each element
# of the list describes a "system"; a system is characterized
# by a parm7 file. If you have multiple conformations and/or
# scans that all use the same parm7 file, then you only have
# 1 system.  In contrast, you may be parametrizing a potential
# that exists in multiple molecules, in which each system
# refers to the same set of training parameters.
#
# A system's dictionary contains several keys.
# parm: the name of an amber parm7 input file.
# crd: the name of a formatted rst7 input file.
# output: the name of a python output file.
# params: a dictionary that describes how to map
# the parameters to the atoms.
# profiles: a list that collects high-level and
# low-level structures used to train the parameters.
# Each profile is a series of structures that
# share an arbitrary, yet common, zero-of-energy.
#
# The each key of the params dictionary is one of
# unique parameters. The value is a list-of-lists.
# Each sublist contains 4 elements: the amber masks
# that select each of the 4 atoms in the torsion
# BY ATOM NAME; e.g., [ "@C1", "@C2", "C3", "@H1" ]
# You do not need to map global parameters; therefore,
# if all parameters were global parameters, then the
# "params" dictionary would be empty, {}. If the same
# unique parameter should be applied to more than one
# quartet, then there would be a sublist for each
# quartet.
#
# The profiles value is a list of dictionaries.
# Each dictionary describes a "scan"; it contains the
# keys: hl, ll, name, and plots.
# The hl value is the name of a scan performed with a
# target model chemistry. These are the energies we
# would like to reproduce.
# The ll value is the scan performed with the force field
# before changing the parameters.
# The name value is a prefix applied to the plots generated
# during the nonlinear optimization procedure.
# The "plots" value is a list of strings used to further
# define filenames.
#
# The following is an example that parametrizes a single
# bespoke potential
#
# {
#    "params": {
#        "param_name_1": {
#           "nprim": 3,
#           "masks": null
#         }
#    },
#    "output": "global.frcmod",
#    "systems": [
#       {
#           "parm": "system1.parm7",
#           "crd":  "system1.rst7",
#           "output": "system1.py",
#           "params": {
#                  "param_name_1": [
#                     [ "@AtomName1", "@AtomName2", "@AtomName3", "@AtomName4" ]
#                  ]
#           },
#           "profiles": [
#               {
#                   "hl": "highlevel_scan.xyz",
#                   "ll": "lowlevel_scan.xyz",
#                   "name": "output_plot_prefix",
#                   "plots": [
#                        "param_name_1"
#                   ]
#               }
#           ]
#       }
# }
#
# The following is an example that optimizes a global parameter.
#
# {
#    "params": {
#        "param_name_1": {
#           "nprim": 3,
#           "masks": [ ["@%nf","@%ce","@%ca","@%ca"] ]
#         }
#    },
#    "output": "global.frcmod",
#    "systems": [
#       {
#           "parm": "system1.parm7",
#           "crd":  "system1.rst7",
#           "output": "system1.py",
#           "params": {},
#           "profiles": [
#               {
#                   "hl": "highlevel_scan.xyz",
#                   "ll": "lowlevel_scan.xyz",
#                   "name": "output_plot_prefix",
#                   "plots": [
#                        "param_name_1"
#                   ]
#               }
#           ]
#       }
# }
#
