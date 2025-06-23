#!/usr/bin/env python3
import sys
import argparse
from parmed import load_file
from parmed.tools.actions import deleteDihedral, addDihedral
from parmed.amber.mask import AmberMask
parser = argparse.ArgumentParser("replace dihedral parameters")
parser.add_argument("--resname",default="LIG",help="Name of the residue, default: LIG",type=str)
parser.add_argument("iparm",help="Input parm7")
parser.add_argument("oparm",help="Output parm7")
args = parser.parse_args()
rname = args.resname
scee = 1.2
scnb = 2.0
if args.iparm == args.oparm:
    raise Exception("The 2 filenames must be different")

p = load_file( args.iparm )

mask=f":{rname}@C2"
res = [i for i in AmberMask(p,mask).Selected()]
if len(res) == 0:
    raise Exception(f"No atoms matching {mask}")
mask=f":{rname}@C3"
res = [i for i in AmberMask(p,mask).Selected()]
if len(res) == 0:
    raise Exception(f"No atoms matching {mask}")
mask=f":{rname}@C4"
res = [i for i in AmberMask(p,mask).Selected()]
if len(res) == 0:
    raise Exception(f"No atoms matching {mask}")
mask=f":{rname}@C7"
res = [i for i in AmberMask(p,mask).Selected()]
if len(res) == 0:
    raise Exception(f"No atoms matching {mask}")
mask=f":{rname}@O7"
res = [i for i in AmberMask(p,mask).Selected()]
if len(res) == 0:
    raise Exception(f"No atoms matching {mask}")
mask=f":{rname}@N8"
res = [i for i in AmberMask(p,mask).Selected()]
if len(res) == 0:
    raise Exception(f"No atoms matching {mask}")
mask=f":{rname}@C9"
res = [i for i in AmberMask(p,mask).Selected()]
if len(res) == 0:
    raise Exception(f"No atoms matching {mask}")
mask=f":{rname}@S10"
res = [i for i in AmberMask(p,mask).Selected()]
if len(res) == 0:
    raise Exception(f"No atoms matching {mask}")
mask=f":{rname}@N13"
res = [i for i in AmberMask(p,mask).Selected()]
if len(res) == 0:
    raise Exception(f"No atoms matching {mask}")


deleteDihedral(p,f":{rname}@S10",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7").execute()
addDihedral(p,f":{rname}@S10",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",-0.6104380159930654,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@S10",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",-2.484919749058815,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@S10",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",-0.6869450510866189,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@N13",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7").execute()
addDihedral(p,f":{rname}@N13",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",1.2690037710703148,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N13",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",-1.7967395818879786,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N13",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",0.32936112937894635,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@C3").execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",0.7787614874785769,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",-5.107378515958192,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",-0.056037412113724545,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@O7").execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@O7",1.9024956376089779,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@O7",-3.3346544708199044,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@O7",-0.9950435977166064,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2").execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",0.32935505210203775,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",-0.31149635696913586,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",-0.07052305978537719,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4").execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",0.32935505210203775,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",-0.31149635696913586,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",-0.07052305978537719,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2").execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",0.4182428395538059,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",-0.2964983828847052,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",0.02918144015727811,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4").execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",0.4182428395538059,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",-0.2964983828847052,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",0.02918144015727811,3,0,scee,scnb).execute()


p.save(args.oparm,overwrite=True)
