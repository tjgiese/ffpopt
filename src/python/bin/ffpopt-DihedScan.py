#!/usr/bin/env python3

if __name__ == "__main__":
    from ffpopt.GeomOpt import FwdRevDihedScan
    from ffpopt.AmberParm import parmed2ase,parmed2graph
    from ffpopt.Options import AddStandardOptions,StandardArgs
    from ffpopt.Constraints import Constraint
    from parmed import load_file
    import ase.io
    import numpy as np
    import argparse
    from pathlib import Path
    from ffpopt.constants import AU_PER_KCAL_PER_MOL
    from ffpopt.constants import AU_PER_ELECTRON_VOLT

    parser = argparse.ArgumentParser \
        ( formatter_class=argparse.RawDescriptionHelpFormatter,
          description="""Perform a relaxed dihedral scan""" )
    
    parser.add_argument \
        ("-d","--dihed",
         help="4 comma-separated 0-based integers defining the torsion",
         required=True,
         type=str)
    
    parser.add_argument \
        ("-D","--delta",
         help="Scan delta (degrees). Default: 10",
         default=10,
         type=int)
    
    parser.add_argument \
        ("-o","--oscan",
         help="Output extxyz file",
         required=True,
         type=str)
    
    parser.add_argument \
        ("--constrain",
         help="comma-separated list of 2,3,or 4 0-based integers. These are extra constraints that are NOT scanned.",
         required=False,
         type=str,
         action='append')

    AddStandardOptions(parser)

    args = parser.parse_args()
    stdargs = StandardArgs(args)

    con = Constraint.from_str(args.dihed,graph=stdargs.graph)
    con.value = None

    extracons = []
    if args.constrain is not None:
        for c in args.constrain:
            extracons.append( Constraint.from_str(c,graph=stdargs.graph) )
    
    if args.delta < 1:
        raise Exception("--delta must be > 0")
    
    sched = np.arange(0,360,args.delta)
    
    mingeom,scan = FwdRevDihedScan(stdargs.atoms,stdargs,con,extracons,sched)
    
    ase.io.write(args.oscan,scan,format="extxyz")

    KCAL_PER_EV = AU_PER_ELECTRON_VOLT() / AU_PER_KCAL_PER_MOL()

    angles=[]
    enes = []
    for s in scan:
        angles.append(s.info["angle"])
        enes.append(s.info["energy"] * KCAL_PER_EV)
    enes = np.array(enes)
    enes = enes - np.amin(enes)
    
    dat = Path(args.oscan).with_suffix(".dat")
    fh = open(dat,"w")
    for a,e in zip(angles,enes):
        fh.write("%12.2f %20.12e\n"%(a,e) )
    fh.close()


    
