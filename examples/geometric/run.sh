#!/bin/bash
set -e
set -u

geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"sander","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="sander" input.xyz cons.10-9-8-6.inp
#geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"mace","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="mace" input.xyz cons.10-9-8-6.inp
#geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"qdpi2","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="qdpi2" input.xyz cons.10-9-8-6.inp
#geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"xtb","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="xtb" input.xyz cons.10-9-8-6.inp
#geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"HF/6-31G*","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="psi4" input.xyz cons.10-9-8-6.inp


