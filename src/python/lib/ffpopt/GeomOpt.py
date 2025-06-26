#!/usr/bin/env python3


def GeomOpt_ASE(atoms,stdargs,constraints=None):
    import os
    import copy
    import subprocess as subp
    import ase.io
    from . constants import AU_PER_ELECTRON_VOLT
    from . Options import argparse2geometric
    from . Constraints import FillConstraints,ApplyConstraints
    from . Constraints import constraints2info,constraints2ase
    
    if True:
        
        from ase.optimize import BFGS
        
        if stdargs.calc is None:
            stdargs.calc = stdargs.MakeCalc()

        cons = copy.deepcopy(constraints)
        myatoms = atoms.copy()
        asecons = None
        if constraints is not None:
            cons = FillConstraints(myatoms,cons)
            myatoms = ApplyConstraints(myatoms,cons)
            asecons = constraints2ase(cons)
            
        del myatoms.constraints
        myatoms.set_constraint( asecons )
        myatoms.calc = stdargs.calc
        myatoms.calc.reset()
        optimizer = BFGS(myatoms)
        optimizer.run(fmax=stdargs.args.ase_opt_tol)

        out = myatoms.copy()
        out.info = {"energy": myatoms.get_potential_energy()}
        constraints2info(out,cons)

    return out



def GeomOpt_GEOMETRIC(atoms,stdargs,constraints=None):
    import os
    import copy
    import subprocess as subp
    import ase.io
    from . constants import AU_PER_ELECTRON_VOLT
    from . Options import argparse2geometric
    from . Constraints import FillConstraints,ApplyConstraints
    from . Constraints import constraints2info,constraints2ase
    from tempfile import mkstemp
    import os
 
    if True:
        #
        # Build the list of geomtric-optimize command line arguments
        #
        parm = stdargs.parm
        try:
            if stdargs.calc is not None:
                #print("stdargs has calc")
                parm = stdargs.calc.parm
                #print("stdargs.calc has parm",parm)
        except:
            pass


        
        fd,tmpxyz = mkstemp(dir=".",prefix="tmp.",suffix=".xyz")
        if not os.isatty(fd):  # Check if fd is still valid
            os.close(fd)

        tmpopt = tmpxyz.replace(".xyz","_optim.xyz")
        tmplog = tmpxyz.replace(".xyz",".log")
        tmpcons = tmpxyz.replace(".xyz",".cons.inp")
        tmpdir = tmpxyz.replace(".xyz",".tmp")

        cmds = argparse2geometric(stdargs.model,parm,stdargs.crd,stdargs.args)
        cmds.append( tmpxyz )
    
        cons = copy.deepcopy(constraints)
        myatoms = atoms.copy()
        if constraints is not None:
            cons = FillConstraints(myatoms,cons)
            myatoms = ApplyConstraints(myatoms,cons)

            cmds.append(tmpcons)
            fh = open(tmpcons,"w")
            fh.write("$set\n")
            for con in cons:
                line = con.to_geometric()
                fh.write("%s\n"%(line))
            fh.close()
        
        ase.io.write(tmpxyz,myatoms)

        print("cmds=",cmds)
    
        result = subp.run(cmds,capture_output=False,text=True)

        if not os.path.exists(tmpopt):
            raise Exception(f"File not found: {tmpopt}")
    
        out = ase.io.read(tmpopt,index='-1')
        keys = [ key for key in out.info ]
        ene = float(keys[-1]) / AU_PER_ELECTRON_VOLT()
        out.info = { "energy": ene }
        constraints2info(out,cons)

        for f in [tmpxyz,tmpopt,tmplog,tmpcons]:
            if os.path.exists(f):
                os.remove(f)

        if os.path.isdir(tmpdir):
            import shutil
            shutil.rmtree(tmpdir)
    
    return out


def GeomOpt_SinglePoint(atoms,stdargs,constraints=None):
    import os
    import copy
    import subprocess as subp
    import ase.io
    from . constants import AU_PER_ELECTRON_VOLT
    from . Options import argparse2geometric
    from . Constraints import FillConstraints,ApplyConstraints
    from . Constraints import constraints2info,constraints2ase

    if True:
        if stdargs.calc is None:
            stdargs.calc = stdargs.MakeCalc()
        out = atoms.copy()
        del out.constraints
        out.calc = stdargs.calc
        out.calc.reset()
        out.info = {"energy": out.get_potential_energy()}
        out.calc = None
        
    return out



def GeomOpt(atoms,stdargs,constraints=None):

    if stdargs.args.no_opt:
        out = GeomOpt_SinglePoint(atoms,stdargs,constraints)
    elif stdargs.args.ase_opt:
        try:
            out = GeomOpt_ASE(atoms,stdargs,constraints)
        except Exception as e:
            import traceback
            print("\n\n\nASE GEOMETRY OPTIMIZATION FAILURE\n")
            print(e)
            traceback.print_exc()
            out = GeomOpt_GEOMETRIC(atoms,stdargs,constraints)
    else:
        out = GeomOpt_GEOMETRIC(atoms,stdargs,constraints)
    return out
            




def ApplyDihedConstraint(atoms,idxs,value,rotmask):
    import copy
    out = atoms.copy()
    out.set_dihedral(idxs[0],idxs[1],idxs[2],idxs[3],value,
                     mask=rotmask)
    return out


    

def DihedScan(atoms,stdargs,con,extracons,sched):
    import numpy as np
    import copy
    
    idxs = copy.deepcopy( con.idxs )

    mingeom = GeomOpt(atoms,stdargs,constraints=extracons)

    sched = np.array(sched,copy=True)
    minang = mingeom.get_dihedral(idxs[0],idxs[1],idxs[2],idxs[3])

    difangs = [ abs((x-minang)%360) for x in sched ]
    istart = np.argmin(difangs)
    nscan = len(sched)
    scan = []

    for i in range(nscan+1):
        icur = (istart + i) % nscan
        if icur == nscan:
            value = sched[0]
        else:
            value = sched[icur]

        if i == 0:
            igeom = mingeom.copy()
        else:
            igeom = opt.copy()

        fang = igeom.get_dihedral(idxs[0],idxs[1],idxs[2],idxs[3])
            
        cons = [ copy.deepcopy(con) ]
        cons[0].value = value

        opt = GeomOpt(igeom,stdargs,constraints=cons + extracons)
        

        print("Finished angle:",value)
        
        opt.info["angle"] = value
        if i == nscan:
            if opt.info["energy"] < scan[0].info["energy"]:
                scan[0] = opt
        else:
            scan.append(opt)

    scan = sorted(scan,key=lambda x: x.info["angle"])

    for s in scan:
        if s.info["energy"] < mingeom.info["energy"]:
            mingeom = s
    
    return mingeom,scan


def FwdRevDihedScan_worker(args):
    return DihedScan(*args)


def FwdRevDihedScan(atoms,stdargs,con,extracons,sched,parallel=False):
    
    import ase.io
    
    revsched = sched[::-1]

    if parallel:
        import multiprocessing

        proclist = [ (atoms,stdargs,con,extracons,sched),
                     (atoms,stdargs,con,extracons,revsched) ]
        
        with multiprocessing.Pool(processes=2) as pool:
            olists = pool.map(FwdRevDihedScan_worker,proclist)
        fmingeom,fwd = olists[0]
        rmingeom,rev = olists[1]

    else:
        fmingeom,fwd = DihedScan(atoms,stdargs,con,extracons,sched)
        rmingeom,rev = DihedScan(atoms,stdargs,con,extracons,revsched)
    
    #ase.io.write("tmp.fscan.extxyz",fwd)
    #ase.io.write("tmp.rscan.extxyz",rev)
    
    scan = []
    
    for a,b in zip(fwd,rev):
        if abs(a.info["angle"] - b.info["angle"]) > 1.e-4:
            raise Exception("Expected fwd and rev scans at the same set of dihedrals, but found %s and %s\n"%(a.info["angle"],b.info["angle"]))
        if a.info["energy"] < b.info["energy"]:
            scan.append(a)
        else:
            scan.append(b)

    if fmingeom.info["energy"] < rmingeom.info["energy"]:
        mingeom = fmingeom
    else:
        mingeom = rmingeom

    return mingeom,scan

