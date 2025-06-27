#!/usr/bin/env python3

#
# The xtb-python package is no longer maintained and the developers
# recommend using tblite instead.
#
# To use the old XTB implementation, install the following packages:
#   conda install -y xtb xtb-python
#
# Then uncomment the following lines:
#
#from xtb.ase.calculator import XTB as OldXTB
#class XTBCalculator(OldXTB):
#     """ASE calculator for XTB with net charge."""
#     def __init__(self, *args, charge=0, **kwargs):
#         self.charge = charge
#         super().__init__(*args, **kwargs)
#     def _create_api_calculator(self):
#         import numpy as np
#         initial_charges = np.zeros(len(self.atoms))
#         initial_charges[0] = self.charge
#         self.atoms.set_initial_charges(initial_charges)
#         return super()._create_api_calculator()
#
#
# To use the new tblite implementation, install the following packages:
#   python3 -m pip install tblite
#
# Then load the calculator:
#  from tblite.ase import TBLite
#  calc = TBLite(method="GFN2-xTB",charge=self.charge,verbosity=-1)
#

from ase.calculators.calculator import Calculator, all_changes
from collections import defaultdict as ddict

def CopyParm( parm ):
    import copy
    try:
        parm.remake_parm()
    except:
        pass
    p = copy.copy( parm )
    p.coordinates = copy.copy( parm.coordinates )
    p.box = copy.copy( parm.box )
    try:
        p.hasbox = copy.copy( parm.hasbox )
    except:
        p.hasbox = False
    return p


class GenCalculator(Calculator):

    implemented_properties = ['energy','forces','free_energy']
    nolabel=True

    def __init__(self,mode="sander",parm="tmp.parm7",crd="tmp.rst7",mol=None,**kwargs):
        from parmed import load_file
        self.parm=parm
        self.crd=crd
        if mol is not None:
            self.mol = mol
        else:
            self.mol = load_file(parm,xyz=crd)
        self.charge = int(round(sum([a.charge for a in self.mol.atoms])))
        self.mode = mode.upper()
        if self.mode == "SANDER":
            self.calc = SanderCalculator(parm=parm,crd=crd,**kwargs)
        elif self.mode == "QDPI2":
            import importlib
            from pathlib import Path
            data_file_name = "pkgdata/qdpi/qdpi-2.0.pb"
            data_path = importlib.resources.files("ffpopt") / data_file_name
            #file_path = importlib.resources.as_file(data_path)
            model = kwargs.get("model",data_path)
            mlp = DPModel(model)
            self.calc = QDpi2Calculator(mlp,self.charge,**kwargs)
        elif self.mode == "XTB":
            #self.calc = XTBCalculator(charge=self.charge,method="GFN2-XTB")
            from tblite.ase import TBLite
            self.calc = TBLite(method="GFN2-xTB",charge=self.charge,verbosity=-1)
        elif self.mode == "MACE":
            from mace.calculators import MACECalculator
            import importlib
            from pathlib import Path
            data_file_name = "pkgdata/mace-off/mace_off23/MACE-OFF23_medium.model"
            data_path = importlib.resources.files("ffpopt") / data_file_name
            #file_path = importlib.resources.as_file(data_path)
            model = kwargs.get("model",data_path)
            #print(model,data_path)
            self.calc = MACECalculator(model_paths=model)
        elif "AIMNET" in self.mode:
            from aimnet.calculators import AIMNet2ASE
            self.calc = AIMNet2ASE(base_calc=self.mode.lower(),charge=self.charge)
        elif "ANI" in self.mode:
            import torchani.models
            if self.mode == "ANI1CCX":
                self.calc = torchani.models.ANI1ccx().ase()
            elif self.mode == "ANI1X":
                self.calc = torchani.models.ANI1x().ase()
            elif self.mode == "ANI2X":
                self.calc = torchani.models.ANI2x().ase()
            else:
                raise Exception(f"Expected ani1x, ani2x, or ani1ccx, but received {self.mode}")
        elif "/" in self.mode:
            from ase.calculators.psi4 import Psi4
            
            import os
            import sys
            cwd = os.getcwd()
            if "PSI_SCRATCH" not in os.environ:
                sys.stderr.write(f"PSI_SCRATCH is unset. Setting to {cwd}\n")
                os.environ["PSI_SCRATCH"] = cwd
            else:
                sdir = os.environ["PSI_SCRATCH"]
                if len(sdir) == 0:
                    sys.stderr.write(f"PSI_SCRATCH='{sdir}' does not refer to a directory. Setting to {cwd}\n")
                    os.environ["PSI_SCRATCH"] = cwd
                else:
                    if not os.path.isdir(sdir):
                        sys.stderr.write(f"PSI_SCRATCH='{sdir}' does not exist. Setting to {cwd}\n")
                        os.environ["PSI_SCRATCH"] = cwd

            
            theory,basis = self.mode.split("/")
            memory = kwargs.get("memory","1gb")
            num_threads = int(kwargs.get("num_threads",4))
            self.calc = Psi4(method=theory,memory=memory,basis=basis,num_threads=num_threads,charge=self.charge,multiplicity=1)
        else:
            raise Exception(f"Unknown mode: {self.mode}")
        Calculator.__init__(self,**kwargs)

    def calculate(self,
                  atoms=None,
                  properties=None,
                  system_changes=all_changes):
        import numpy as np
        import ase
        if properties is None:
            properties = self.implemented_properties
        Calculator.calculate(self, atoms, properties, system_changes)
        atoms.calc =  self.calc
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        #print(f"Energy = {energy} eV")
        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces






class DPModel(object):
    
    def __init__(self,fname):
        
        from deepmd.infer import DeepPot
        #from deepmd.env import reset_default_tf_session_config

        #try:
        self.dp = DeepPot(fname)
        #except:
        #    reset_default_tf_session_config(True)
        #    self.dp = DeepPot(fname)

        self.cell   = None
        self.rcut   = self.dp.get_rcut()
        self.ntypes = self.dp.get_ntypes()
        self.tmap   = self.dp.get_type_map()

    def GetTypeIdxFromSymbol(self,ele):
        idx = None
        if ele in self.tmap:
            idx = self.tmap.index(ele)
        return idx

    def GetTypeIdxs(self,eles):
        return [ self.GetTypeIdxFromSymbol(ele) for ele in eles ]

    def CalcEne(self,eles,crds):
        import numpy as np
        
        #from dpdata.unit import EnergyConversion
        #from dpdata.unit import ForceConversion
        #from dpdata.unit import LengthConversion

        # Bohr/angstrom
        #length_convert = LengthConversion("angstrom", "bohr").value()
        # hartree/eV
        #energy_convert = EnergyConversion("eV", "hartree").value()
        # hartree/bohr / eV/angstrom
        #force_convert = ForceConversion("eV/angstrom", "hartree/bohr").value()

        energy_convert = 1
        force_convert = 1
        
        coord = np.array(crds).reshape([1, -1])
        atype = self.GetTypeIdxs(eles)
        e, f, v = self.dp.eval(coord, self.cell, atype)
        f = f[0] * force_convert
        e = e[0][0] * energy_convert

        # print("$coord")
        # for i in range(len(eles)):
        #     print("%20.14f %20.14f %20.14f %s"%(
        #         crds[i,0] * length_convert,
        #         crds[i,1] * length_convert,
        #         crds[i,2] * length_convert,
        #         eles[i]))
        # print("$end")
        
        return e,f


    

    


class QDpi2Calculator(Calculator):

    implemented_properties = ['energy','forces','free_energy']
    nolabel=True
    
    def __init__(self,dpmodel,charge,**kwargs):
        from tblite.ase import TBLite
        self.dpmodel = dpmodel
        self.charge = charge
        #self.xtbcalc = XTBCalculator(charge=charge,method="GFN2-xTB")
        self.xtbcalc = TBLite(method="GFN2-xTB",charge=self.charge)
        Calculator.__init__(self,**kwargs)
        
    def calculate(self,
                  atoms=None,
                  properties=None,
                  system_changes=all_changes):
        import numpy as np
        import ase
        if properties is None:
            properties = self.implemented_properties
        Calculator.calculate(self, atoms, properties, system_changes)
        natoms = len(self.atoms)
        positions = self.atoms.positions
        forces = np.zeros((natoms, 3))
        energy = 0

        eles = self.atoms.get_chemical_symbols()
        crds = self.atoms.get_positions()
        
        e,f = self.dpmodel.CalcEne(eles,crds)
        
        atlist = "".join( ["%s1"%(ele) for ele in eles ] )
        atoms = ase.Atoms(atlist,positions=crds)
        atoms.calc =  self.xtbcalc
        e2 = atoms.get_potential_energy()
        f2 = atoms.get_forces()

        energy = e + e2
        forces = f + f2

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces


class SanderCalculator(Calculator):

    implemented_properties = ['energy','forces','free_energy']
    nolabel=True

    def __init__(self,parm="tmp.parm7",crd="tmp.rst7",mol=None,**kwargs):
        from ase.calculators.amber import SANDER
        import sander
        from parmed import load_file
        self.parm = parm
        if mol is not None:
            self.mol = CopyParm(mol)
        else:
            self.mol = load_file(parm,xyz=crd)
        self.mm_options = sander.gas_input()
        self.mm_options.cut=99.
        self.mm_options.ntc=1
        self.mm_options.ntf=1
        self.calc = SANDER(top=self.parm,crd=self.mol,mm_options=self.mm_options)
        Calculator.__init__(self,**kwargs)

    def calculate(self,
                  atoms=None,
                  properties=None,
                  system_changes=all_changes):
        import numpy as np
        import ase
        if properties is None:
            properties = self.implemented_properties
        Calculator.calculate(self, atoms, properties, system_changes)
        atoms.calc =  self.calc
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces




