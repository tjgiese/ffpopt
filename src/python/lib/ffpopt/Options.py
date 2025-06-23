#!/usr/bin/env python3

def AddGeomOptOptions(parser):
   
   parser.add_argument \
      ("--no-opt",
       help="Skip all geometry optimizations and simply perform an energy evaluation",
       action='store_true')
   
   parser.add_argument \
      ("--ase-opt",
       help="If present, then use the BFGS optimizer in ASE rather than geometric-optimize",
       action='store_true')

   parser.add_argument \
      ("--ase-opt-tol",
       help="The ASE geometry optimization tolerance. Default: 0.01 eV/A",
       default=0.01,
       type=float)
   
   parser.add_argument \
       ("--geometric-maxiter",
        help="Maximum number of optimization steps. Default: 500",
        default=500,
        type=int)
   
   parser.add_argument \
       ("--geometric-coordsys",
        help="Coordinate system. Default: tric",
        default='tric',
        type=str)
   
   parser.add_argument \
       ("--geometric-converge",
        help="Optimization tolerance(s). Default: 'set GAU_TIGHT'",
        default='set GAU_TIGHT',
        type=str)

   parser.add_argument \
       ("--geometric-enforce",
        help="Constraint enforcement tolerance. Default: 0.1",
        default=0.1,
        type=float)


def AddModelOptions(parser):
   
   parser.add_argument \
       ("-m","--model",
        help="Energy calculator: sander, xtb, qdpi2, mace, aimnet2, aimnet2_wb97m (aimnet2 is an alias for aimnet2_wb97m), aimnet2_b973c, aimnet2_qr, ani1x, ani2x, ani1ccx, theory/basis (psi4). Default: sander",
        default="sander",
        type=str)

   parser.add_argument \
       ("--psi4-memory",
        help="The available memory. This is sent to the psi4 calculator. Default: '1gb'",
        default='1gb',
        type=str)

   parser.add_argument \
       ("--psi4-num-threads",
        help="The number of threads sent to the psi4 calculator. Default: 4",
        default=4,
        type=int)


def AddSanderOptions(parser):
   parser.add_argument \
       ("-p","--parm",
        help="Input Amber parm7",
        required=True,
        type=str)

   parser.add_argument \
       ("-c","--crd",
        help="Input Amber rst7",
        required=True,
        type=str)



def AddStandardOptions(parser):
   AddSanderOptions(parser)
   AddModelOptions(parser)
   AddGeomOptOptions(parser)
   


def GetStandardOptions(args):

   geometric_kwargs = {}
   try:
      geometric_kwargs = { "coordsys": str(args.geometric_coordsys),
                           "maxiter": str(args.geometric_maxiter),
                           "converge": str(args.geometric_converge),
                           "enforce": str(args.geometric_enforce) }
   except:
      pass

   
   psi4_kwargs = {}
   try:
      psi4_kwargs = { "memory": str(args.psi4_memory),
                      "num_threads": str(args.psi4_num_threads) }
   except:
      pass

   
   extra_args = { "geometric": geometric_kwargs,
                  "psi4": psi4_kwargs }

   return extra_args


class StandardArgs(object):

   @classmethod
   def from_manual_parm(cls,parm,crd,args):
      import copy
      fargs = copy.deepcopy(args)
      fargs.parm = parm
      fargs.crd = crd
      return cls(fargs)

   
   def __init__(self,args):
      
      from parmed import load_file
      from . AmberParm import parmed2ase, parmed2graph
      
      self.args = args
      self.parm = args.parm
      self.crd = args.crd
      self.mol = load_file(self.parm,xyz=self.crd)
      self.atoms, self.charge = parmed2ase(self.mol)
      self.graph = parmed2graph(self.mol)
      self.model = args.model
      self.calc = None

      
   def MakeCalc(self,parm=None,crd=None):
      
      from . ase.calculator import GenCalculator
      
      if parm is None:
         parm = self.args.parm
      if crd is None:
         crd = self.args.crd

      memory = getattr(self.args,"psi4_memory","1gb")
      num_threads = getattr(self.args,"psi4_num_threads",1)
         
      kwargs = { "memory": memory,
                 "num_threads": num_threads }

      try:
         import sander
         sander.cleanup()
      except:
         pass
         
      return GenCalculator(mode=self.args.model,
                           parm=parm,
                           crd=crd,
                           mol=self.mol,
                           **kwargs)
                           
         


def ModelIsPsi4(model):
   m = model.upper()
   ispsi4 = False
   if ".pb" in m:
      pass
   elif ".model" in m:
      pass
   elif "XTB" in m:
      pass
   elif "QDPI2" in m:
      pass
   elif "MACE" in m:
      pass
   elif "AIMNET" in m:
      pass
   elif "/" in m:
      ispsi4 = True
   return ispsi4


def argparse2geometric(model,parm,crd,args):
    
   kwargs = GetStandardOptions(args)
   geok = kwargs["geometric"]
   psik = kwargs["psi4"]

   geoopts = []
   for key in geok:
      val = geok[key]
      geoopts.append( "--%s"%(key) )
      if key == "converge":
         geoopts.extend( geok[key].split() )
      else:
         geoopts.append( geok[key] )

   asek = { "mode": model, "parm": parm, "crd": crd }
   if ModelIsPsi4(model):
      for key in psik:
         asek[key] = str(psik[key])
   
            
   asestr = ",".join( [ '"%s": "%s"'%(key,asek[key]) for key in asek ])
   asestr = "{%s}"%(asestr)
            
   cmds = ["geometric-optimize",
           "--engine", "ase",
           "--ase-class", "ffpopt.ase.GenCalculator",
           "--ase-kwargs", asestr ] + geoopts
    
   return cmds
