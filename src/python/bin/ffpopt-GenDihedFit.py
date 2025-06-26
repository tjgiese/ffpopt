#!/usr/bin/env python3


            
if __name__ == "__main__":
    import argparse
    import numpy as np
    from ffpopt.Options import AddGeomOptOptions
    from ffpopt.Dihedrals import FitInputType
    from ffpopt.Dihedrals import NonlinearSolve

    parser = argparse.ArgumentParser \
        ( formatter_class=argparse.RawDescriptionHelpFormatter,
          description="""Perform a general, multidihedral fit""" )

    parser.add_argument\
        ("inp",
         help="Json input file")

    parser.add_argument\
        ("--stride",
         default=1,
         type=int,
         help="Stride when reading structures. Default: 1")
    
    parser.add_argument \
        ("--nlmaxiter",
         help="Maximum number of nonlinear optimization steps. Default: 300",
         default=300,
         type=int)

    parser.add_argument \
        ("--nlrhobeg",
         help="Initial parameter displacements. Default: 0.25 kcal/mol",
         default=0.25,
         type=float)
    
    
    parser.add_argument \
        ("--nltol",
         help="Tolerance on the parameter optimization. Default: 0.01",
         default=0.01,
         type=float)
    

    AddGeomOptOptions(parser)
    args = parser.parse_args()
    args.model="sander"
    finp = FitInputType.from_file(args,args.inp)
    NonlinearSolve(args,finp)
    finp.write_output()
