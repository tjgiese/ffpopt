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
addDihedral(p,f":{rname}@S10",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",-0.7803437597780679,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@S10",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",-2.6680734399571966,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@S10",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",-0.3228598091433068,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@N13",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7").execute()
addDihedral(p,f":{rname}@N13",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",1.0442375322237256,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N13",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",-1.7605963110451317,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N13",f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",0.6433317105272123,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@C3").execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",0.04314905450282817,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",-3.675147650272786,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",0.7617627641030837,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@O7").execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@O7",1.7195490929679744,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@O7",-4.104654145499766,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@C9",f":{rname}@N8",f":{rname}@C7",f":{rname}@O7",-1.5947036922460525,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2").execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",0.7136919035917111,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",-0.5311181223138036,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",-0.07634995027444282,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4").execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",0.7136919035917111,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",-0.5311181223138036,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@O7",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",-0.07634995027444282,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2").execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",0.1515228704325982,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",0.2963098289722953,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C2",0.874378619794888,3,0,scee,scnb).execute()


deleteDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4").execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",0.1515228704325982,1,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",0.2963098289722953,2,0,scee,scnb).execute()
addDihedral(p,f":{rname}@N8",f":{rname}@C7",f":{rname}@C3",f":{rname}@C4",0.874378619794888,3,0,scee,scnb).execute()


p.save(args.oparm,overwrite=True)
