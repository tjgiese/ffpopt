#!/usr/bin/env python3

if __name__ == "__main__":
    from ffpopt.GeomOpt import GeomOpt
    from ffpopt.AmberParm import parmed2ase, parmed2graph
    from ffpopt.Options import AddStandardOptions,StandardArgs
    from ffpopt.Constraints import Constraint
    from ffpopt.Reader import ReadGeomsFromXYZ
    from parmed import load_file
    import ase.io
    import numpy as np
    import argparse
    
    parser = argparse.ArgumentParser \
        ( formatter_class=argparse.RawDescriptionHelpFormatter,
          description="""Perform a geometry optimization""" )
    
    parser.add_argument \
        ("--constrain",
         help="comma-separated list of 2,3,or 4 0-based integers. The list can be appended with =value to specify a value. If not given, then the input coordinates are used to assign the value. This can be used multiple times.",
         required=False,
         type=str,
         action='append')
    
    parser.add_argument \
        ("-o","--oscan",
         help="Output extxyz file",
         required=True,
         type=str)

    parser.add_argument \
        ("-i","--iscan",
         help="Input extxyz file",
         required=False,
         type=str)

    parser.add_argument \
        ("-I","--ignore",
         help="If present, then ignore the constraint information within iscan",
         action='store_true')


    
    
    AddStandardOptions(parser)

    args = parser.parse_args()
    stdargs = StandardArgs(args)

    cons = None
    if args.constrain is not None:
        if len(args.constrain) > 0:
            cons = [ Constraint.from_str(c) for c in args.constrain ]
            for c in cons:
                c.SetMask(stdargs.graph)

    if args.iscan is not None:
        geoms,mcons = ReadGeomsFromXYZ(stdargs,args.iscan)
        if args.ignore:
            mcons = [None]*(len(mcons))
        for atoms,econ in zip(geoms,mcons):
            if econ is not None and cons is not None:
                raise Exception("Input structure has constraints, "+
                                "but constraints were also given on "+
                                "the command line. Use --ignore")
        oxyz = []
        for atoms,econ in zip(geoms,mcons):
            if econ is None:
                econ = cons
            s = GeomOpt(atoms,stdargs,constraints=econ)
            oxyz.append(s)
        ase.io.write(args.oscan,oxyz,format="extxyz")
        
    else:
        s = GeomOpt(stdargs.atoms,stdargs,constraints=cons)
        ase.io.write(args.oscan,s,format="extxyz")



