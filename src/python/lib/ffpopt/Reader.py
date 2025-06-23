#!/usr/bin/env python3

def ReadGeomsFromXYZ(stdargs,fname):
    import ase.io
    from . Constraints import info2constraints, FillConstraints
    
    
    geoms = ase.io.read(fname,index=":")
    cons = []
    for atoms in geoms:
        try:
            atoms.info["energy"] = atoms.get_potential_energy()
        except:
            pass
        mycons=None
        if "constraints" in atoms.info:
            mycons = info2constraints(atoms,stdargs.graph)
            mycons = FillConstraints(atoms,mycons)
        cons.append(mycons)

    return geoms,cons
