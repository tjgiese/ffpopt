#!/usr/bin/env python3


def MakeUniqueDihedralScans(stdargs,scans):
    from collections import defaultdict as ddict

    nscan = len(scans)
    scanmap = ddict(list)
    for i in range(nscan):
        key = scans[i][0].key
        scanmap[key].append(i)

            
    ukeys = [key for key in scanmap]
    #for i in range(nscan):
    #    k = ukeys.index( scans[i][0].key )
    #    scans[i].pidx = k
    
    dprofs = []
    for key in ukeys:
        i = scanmap[key][0]
        con = scans[i][0].con
        profs = [ scans[j] for j in scanmap[key] ]
        iscans = [ j for j in scanmap[key] ]
        dprofs.append( DihedralProfile(stdargs,con,profs,iscans) )

    return dprofs



class DihedralProfile(object):
    def __init__(self,stdargs,con,profs,iscans):
        self.con = con
        self.profs = profs
        self.iscans = iscans
        self.dfcn = IsolatedLinearSolve\
            (stdargs,self.profs[0],self.con,
             stdargs.args.nprim,self.iscans[0])


        
    
class ScanPair(object):
    def __init__(self,stdargs,llgeom,hlgeom):
        from ffpopt.Constraints import info2constraints
        from ffpopt.Constraints import FillConstraints
        import copy

        
        self.llgeom = llgeom.copy()
        self.llcons = info2constraints( self.llgeom, stdargs.graph )
        self.llcons = FillConstraints(self.llgeom,self.llcons,force=True)
        
        self.hlgeom = hlgeom.copy()
        self.hlcons = info2constraints( self.hlgeom, stdargs.graph )
        self.hlcons = FillConstraints(self.hlgeom,self.hlcons,force=True)

        if len(self.llcons) != 1:
            raise Exception("Expected 1 LL constraint")
        
        if len(self.hlcons) != 1:
            raise Exception("Expected 1 HL constraint")

        if len(self.llcons[0].idxs) != 4:
            raise Exception("Expected LL dihedral constraint")
        
        if len(self.hlcons[0].idxs) != 4:
            raise Exception("Expected HL dihedral constraint")

        dang = abs(self.llcons[0].value - self.hlcons[0].value)
        dang = min(dang,360-dang)

        if dang > 0.0001:
            raise Exception("LL and HL constraints differ: "
                            +f"{self.llcons[0].value} vs {self.hlcons[0].value}")
        self.con = copy.deepcopy(self.llcons[0])
        tmp = copy.deepcopy(self.con)
        tmp.value = None
        self.key = str(tmp)
        if stdargs.args.bytype:
            self.key = tuple( [stdargs.mol.atoms[i].type for i in tmp.idxs ] )


def EnergyScansWithoutDihedrals(stdargs,list_of_scans,cons):
    from ffpopt.AmberParm import CopyParm
    from ffpopt.Dihedrals import DeleteDihedrals
    from ffpopt.Dihedrals import GetMultiDihedFcnFromIdxs
    
    p = CopyParm(stdargs.mol)
    
    DeleteDihedrals(p,[ x.idxs for x in cons ])

    p.save("tmp.parm7",overwrite=True)
    calc = stdargs.MakeCalc(parm="tmp.parm7")
    
    list_of_enes = []
    for scan in list_of_scans:
        enes = []
        for geom in scan:
            t = geom.copy()
            del t.constraints
            t.calc = calc
            t.calc.reset()
            e = t.get_potential_energy()
            enes.append(e)
        list_of_enes.append(enes)

    return list_of_enes





def IsolatedLinearSolve(stdargs,scanpairs,con,nprim,CNT):
    from ffpopt.Dihedrals import GetDihedClasses
    from ffpopt.Constraints import FillConstraints
    import numpy as np
    import copy
    from ffpopt.constants import AU_PER_KCAL_PER_MOL
    from ffpopt.constants import AU_PER_ELECTRON_VOLT

    KCAL_PER_EV = AU_PER_ELECTRON_VOLT() / AU_PER_KCAL_PER_MOL()
    
    dfcns = GetDihedClasses(idxs=con.idxs)[nprim]
    llscan = [ x.llgeom for x in scanpairs ]
    list_of_scans = [ llscan ]
    cons = [ con ]
    list_of_enes = EnergyScansWithoutDihedrals(stdargs,list_of_scans,cons)
    #print(list_of_enes[0])
    llenes = np.array(list_of_enes[0])
    hlenes = np.array([ x.hlgeom.info["energy"] for x in scanpairs ])

    llenes *= KCAL_PER_EV
    hlenes *= KCAL_PER_EV

    hlenes -= np.amin(hlenes)
    llenes -= np.amin(llenes)
    
    y = hlenes-llenes

    angs = []
    for igeom,llgeom in enumerate(llscan):
        o = FillConstraints(llgeom,cons,force=True)
        v = o[0].value
        if igeom < 5 and abs(v-360) < 0.00001:
            v = 0
        angs.append( v )
    angs = np.array(angs)
    
    npts = len(y)
    npar = nprim + 1

    bestdfcn = None
    bestchisq = 1.e+30
    bestvalues = []
    for ifcn,dfcn in enumerate(dfcns):
        A = np.zeros( (npts,npar) )
        for iprim,prim in enumerate(dfcn.prims):
            A[:,iprim] = prim.CptEterm(angs)
        A[:,nprim] = 1
        #AtA = A.T @ A
        #AtAinv = np.linalg.inv(AtA)
        #x = AtAinv @ A.T @ y
        x = np.linalg.pinv(A) @ y
        const = x[-1]
        dfcn.SetFCs( x[:-1] )
        v = dfcn.CptEne(angs) + const
        d = hlenes - (llenes + v)
        chisq = np.dot(d,d)
        if chisq < bestchisq:
            bestchisq = chisq
            bestdfcn = copy.deepcopy(dfcn)
            bestvalues = v

    fh = open(f"iso.{CNT}.dat","w")
    for i in range(npts):
        fh.write("%12.3f %20.10e %20.10e %20.10e\n"%\
                 ( angs[i], hlenes[i], llenes[i],
                   llenes[i]+bestvalues[i] ) )
    fh.close()
    
    return bestdfcn
    
    


# def CoupledLinearSolve(stdargs,loscans,cons,dfcns,CNT):
#     from ffpopt.AmberParm import GetDihedClasses
#     from ffpopt.Constraints import FillConstraints
#     import numpy as np
#     import copy
#     from ffpopt.constants import AU_PER_KCAL_PER_MOL
#     from ffpopt.constants import AU_PER_ELECTRON_VOLT

#     KCAL_PER_EV = AU_PER_ELECTRON_VOLT() / AU_PER_KCAL_PER_MOL()

#     llscans = [ [x.llgeom for x in scan ] for scan in loscans ]
    
#     lollenes = EnergyScansWithoutDihedrals(stdargs,llscans,cons)
#     lonpts = []
#     for i in range(len(lollenes)):
#         lollenes[i] = np.array(lollenes[i])
#         lollenes[i] *= KCAL_PER_EV
#         lollenes[i] -= np.amin(lollenes[i])
#         lonpts.append( len(lollenes[i]) )

#     lohlenes = [
#         np.array([ x.hlgeom.info["energy"] for x in scan ])
#         for scan in loscans ]
#     for i in range(len(lohlenes)):
#         lohlenes[i] *= KCAL_PER_EV
#         lohlenes[i] -= np.amin(lohlenes[i])

#     flatll = np.array([ x for scan in lollenes for x in scan ])
#     flathl = np.array([ x for scan in lohlenes for x in scan ])
    
#     y = flathl-flatll
#     tnpts = sum(lonpts)
#     nscan = len(lonpts)



#     loangs = []
#     for iscan in range(nscan):
#         angs = []
#         for igeom,llgeom in enumerate(llscans[iscan]):
#             o = FillConstraints(llgeom,cons)
#             for i in range(len(o)):
#                 v = o[i].value
#                 if igeom < 5 and abs(v-360) < 0.00001:
#                     v = 0
#                 o[i] = v
#             loangs.append(o)
#     loangs = np.array(loangs)
        

    
#     #
#     # List of number of primitive for each dihedral/scan
#     #
    
#     lonprims = [ len(x.prims) for x in dfcns ]
#     npar = sum( lonprims ) + nscan

#     #
#     # Integer offset for each dihedral/scan
#     # The integers are offsets to the force constants
#     # within the flattened parameter array
#     #
    
#     o = 0
#     lopoff = [0]
#     lopidxs = [ [ x for x in range(lonprims[0]) ] ]
#     for iscan in range(1,nscan):
#         o += lonprims[iscan]
#         lopoff.append(o)
#         lopidxs.append( [ x for x in range(o,o+lonprims[iscan]) ] )

#     #
#     # A list of integers for each scan
#     # The integers are the data points within the scan
#     #
    
#     loidxs = []
#     o = 0
#     for iscan in range(nscan):
#         idxs = [x for x in range(o,o+lonpts[iscan])]
#         o += lonpts[iscan]
#         loidxs.append(idxs)

#     A = np.zeros( (tnpts,npar) )

#     #
#     # Loop over data points
#     #
#     for jscan in range(nscan):
#         angs = loangs[jscan]
#         jidxs = loidxs[jscan]
#         #
#         # Loop over force constants
#         #
#         for iscan in range(nscan):
#             o = lopoff[iscan]
#             for iprim,prim in enumerate(dfcns[iscan].prims):
#                 A[jidxs,o+iprim] = prim.CptEterm(loangs[jidxs,iscan])
#         A[jidxs,npar-nscan+jscan] = 1

#     #for ipt in range(tnpts):
#     #    print(" ".join(["%10.2e"%(x) for x in A[ipt,:]]))
        
#     Apinv = np.linalg.pinv(A)
#     x = Apinv @ y
        
#     consts = [ x[npar-nscan+iscan] for iscan in range(nscan) ]
#     #print("x=",x)
#     for iscan in range(nscan):
#         #print("x[iscan]=",x[lopidxs[iscan]])
#         dfcns[iscan].SetFCs( x[lopidxs[iscan]] )
#     #print("consts=",consts)

#     vs = [ np.zeros( (lonpts[iscan],) ) for iscan in range(nscan) ]
    
#     for jscan in range(nscan):
#         for iscan,dfcn in enumerate(dfcns):
#             vs[jscan] += dfcn.CptEne(loangs[loidxs[jscan],iscan])
#         vs[jscan] += consts[jscan]
            
#     chisqs = []
#     for iscan in range(nscan):
#         d = lohlenes[iscan] - (lollenes[iscan] + vs[iscan])
#         chisq = np.dot(d,d)
#         chisqs.append(chisq)
#     chisq = np.sum(chisqs)
    

#     for iscan in range(nscan):
        
#         fh = open(f"mfit.{CNT}.{iscan}.dat","w")
#         angs = loangs[loidxs[iscan],iscan]
#         hlenes = lohlenes[iscan]
#         llenes = lollenes[iscan]
#         cenes = llenes + vs[iscan]
#         for i in range(len(angs)):
#             fh.write("%12.3f %20.10e %20.10e %20.10e\n"%\
#                      ( angs[i], hlenes[i], llenes[i],
#                        cenes[i] ) )
#         fh.close()

#     return chisq,dfcns



class NonlinearObjective(object):
    def __init__(self,stdargs,udscans):
        self.it = 1
        self.stdargs = stdargs
        self.udscans = udscans


        nprims = [ len(udscan.dfcn.prims) for udscan in self.udscans ]
        oprims = [0]
        for iprim in range(1,len(nprims)):
            oprims.append( oprims[-1] + nprims[iprim-1] )
        self.oprims = oprims
        self.parammask = [ [x for x in range(o,o+n)]
                           for o,n in zip(self.oprims,nprims) ]
        self.nparam = sum(nprims)
        self.nscan = len(self.udscans)

        
        
    def GetParams(self):
        import numpy as np
        x = []
        for udscan in self.udscans:
            for prim in udscan.dfcn.prims:
                x.append(prim.fc)
        return np.array(x)


    
    def UpdateFCs(self,x):
        for i in range(self.nscan):
            self.udscans[i].dfcn.SetFCs( x[self.parammask[i]] )


            
    def UpdateFCsAndMakeNewParam(self,params):
        from parmed import DihedralType
        from ffpopt.AmberParm import CopyParm
        from ffpopt.Dihedrals import ChangeDihedrals

        self.UpdateFCs(params)
        p = CopyParm(self.stdargs.mol)
        
        scee = 1.2
        scnb = 2.0
        for iscan in range(self.nscan):
            xs = []
            dfcn = self.udscans[iscan].dfcn
            for prim in dfcn.prims:
                per = prim.per
                ph = prim.phase
                fc = prim.fc
                x = DihedralType(fc, per, ph, scee, scnb)
                xs.append(x)
            ChangeDihedrals(p,dfcn.idxs,xs,bytype=self.stdargs.args.bytype)
        return p

    
            
    def calculate(self,params):
        import numpy as np
        from ffpopt.GeomOpt import GeomOpt
        from ffpopt.constants import AU_PER_KCAL_PER_MOL
        from ffpopt.constants import AU_PER_ELECTRON_VOLT

        KCAL_PER_EV = AU_PER_ELECTRON_VOLT() / AU_PER_KCAL_PER_MOL()
        
        p = self.UpdateFCsAndMakeNewParam(params)
        p.save("tmp.parm7",overwrite=True)
        self.stdargs.calc = self.stdargs.MakeCalc(parm="tmp.parm7")
        x = self.GetParams()
        
        chisq = 0
        iscan = -1
        for udscan in self.udscans:
            for scan in udscan.profs:
                iscan += 1
                llene = []
                hlene = []
                angs = []
                for igeom,pairgeom in enumerate(scan):
                    atoms = pairgeom.llgeom.copy()
                    cons = pairgeom.llcons
                    v = cons[0].value
                    if igeom < len(scan)//2:
                        if abs(v-360) < 0.0001:
                            v = 0
                    angs.append( v )
                    atoms = GeomOpt(atoms,self.stdargs,constraints=cons)
                    hlene.append( pairgeom.hlgeom.info["energy"] * KCAL_PER_EV )
                    llene.append( atoms.info["energy"] * KCAL_PER_EV )
                    if not self.stdargs.args.no_overwrite_guess:
                        pairgeom.llgeom.set_positions( atoms.get_positions() )
                    
                llene  = np.array(llene)
                llene -= np.amin(llene)
                hlene  = np.array(hlene)
                hlene -= np.amin(hlene)
                d = hlene-llene
                llene += np.mean(d)
                d = hlene-llene
                mychisq = np.dot(d,d)
                chisq += mychisq

                hlene -= np.amin(hlene)
                llene -= np.amin(llene)
            
                fh = open(f"nfit.{iscan:02d}.{self.it:05d}.dat","w")
                fh.write("# %25.14f : %s\n"%(mychisq, " ".join(["%.3f"%(w) for w in x])))
                for i in range(len(angs)):
                    fh.write("%20.10e %20.10e %20.10e\n"%\
                             ( angs[i], hlene[i], llene[i] ) )
                fh.close()

        self.it += 1
        return chisq

    
                                    
def ObjFcn(x,self):
    return self.calculate(x)



def NonlinearSolve(stdargs,udscans):
    from scipy.optimize import minimize
    import numpy as np
    import copy

    objfcn = NonlinearObjective(stdargs,udscans)
    n = objfcn.nparam
    x = objfcn.GetParams()
    xlo = x[:] - 2.
    xhi = x[:] + 5.
    #print("n=",n)
    #print("x=",x)
    #print("xlo=",xlo)
    #print("xhi=",xhi)
    bounds = [ (lo,hi) for lo,hi in zip(xlo,xhi) ]
    
    res = minimize( ObjFcn, x, args=(objfcn,),
                    method='COBYLA',
                    bounds=bounds,
                    options={ "rhobeg": stdargs.args.nlrhobeg,
                              "tol": stdargs.args.nltol,
                              "maxiter": stdargs.args.nlmaxiter,
                              "disp": True })

    print(res)

    objfcn.UpdateFCs(res.x)
    #return objfcn.dfcns



def WriteParmedScript(fname,p,dfcns,bytype):
    from collections import defaultdict as ddict
    from ffpopt.Dihedrals import FindDihedrals
    
    aidxs = [ idx for dfcn in dfcns for idx in dfcn.idxs ]
    aidxs = list(set(aidxs))
    resname = p.atoms[aidxs[0]].residue.name
    
    fh = open(fname,"w")
    fh.write("#!/usr/bin/env python3\n")
    fh.write("import sys\n")
    fh.write("import argparse\n")
    fh.write("from parmed import load_file\n")
    fh.write("from parmed.tools.actions import deleteDihedral, addDihedral\n")
    fh.write("from parmed.amber.mask import AmberMask\n")

    fh.write("parser = argparse.ArgumentParser(\"replace dihedral parameters\")\n")
    fh.write(f"parser.add_argument(\"--resname\",default=\"{resname}\",help=\"Name of the residue, default: {resname}\",type=str)\n")
    fh.write("parser.add_argument(\"iparm\",help=\"Input parm7\")\n")
    fh.write("parser.add_argument(\"oparm\",help=\"Output parm7\")\n")
    fh.write("args = parser.parse_args()\n")
    fh.write("rname = args.resname\n")

    fh.write("scee = 1.2\n")
    fh.write("scnb = 2.0\n")

    
    fh.write("if args.iparm == args.oparm:\n")
    fh.write("    raise Exception(\"The 2 filenames must be different\")\n\n")
    
    fh.write("p = load_file( args.iparm )\n\n")
    

    if not bytype:
        for aidx in aidxs:
            res = p.atoms[aidx].residue.name
            name = p.atoms[aidx].name
            mask = ":%s@%s"%("{rname}",name)
            fh.write(f"mask=f\"{mask}\"\n")
            fh.write("res = [i for i in AmberMask(p,mask).Selected()]\n")
            fh.write("if len(res) == 0:\n")
            fh.write("    raise Exception(f\"No atoms matching {mask}\")\n")
            

    fh.write("\n\n")
    for dfcn in dfcns:
        if bytype:
            idxs = dfcn.idxs
            
            ftypes = tuple([p.atoms[idx].type for idx in idxs])
            rtypes = tuple(list(ftypes)[::-1])
            allidxs = []
            for x in p.dihedrals:
                kidxs = (x.atom1.idx,x.atom2.idx,x.atom3.idx,x.atom4.idx)
                tidxs = (x.atom1.type,x.atom2.type,x.atom3.type,x.atom4.type)
                if tidxs == ftypes:
                    allidxs.append(kidxs)
                elif tidxs == rtypes:
                    allidxs.append( tuple(list(kidxs)[::-1]) )
            allmasks = [ [ ":%s@%s"%("{rname}",p.atoms[idx[i]].name) for i in range(4) ]
                         for idx in allidxs ]
        else:
            allmasks = [ [ ":%s@%s"%("{rname}",p.atoms[idx].name)
                           for idx in dfcn.idxs ] ]

        for masks in allmasks:
            mstr = ",".join(["f\"%s\""%(mask) for mask in masks])
            fh.write(f"deleteDihedral(p,{mstr}).execute()\n")
            for prim in dfcn.prims:
                fh.write(f"addDihedral(p,{mstr},{prim.fc},{prim.per},{prim.phase},scee,scnb).execute()\n")
            fh.write("\n\n")

    fh.write("p.save(args.oparm,overwrite=True)\n")
    fh.close()
    
        

if __name__ == "__main__":
    from ffpopt.GeomOpt import GeomOpt
    from ffpopt.AmberParm import parmed2ase, parmed2graph
    from ffpopt.Options import AddStandardOptions,StandardArgs
    from ffpopt.Constraints import Constraint
    from parmed import load_file
    import ase.io
    import numpy as np
    import argparse
    
    parser = argparse.ArgumentParser \
        ( formatter_class=argparse.RawDescriptionHelpFormatter,
          description="""Perform a multidihedral fit""" )
    
    parser.add_argument \
        ("--llscan",
         help="Extended XYZ file of the low-level scan produced by ffpopt-DihedScan.py. Use this option multiple times for multiple dihedrals.  This option must be used the same number of times and in the same order as --hlscan",
         required=True,
         type=str,
         action='append')

    parser.add_argument \
        ("--hlscan",
         help="Extended XYZ file of the high-level scan produced by ffpopt-DihedScan.py. Use this option multiple times for multiple dihedrals.  This option must be used the same number of times and in the same order as --llscan",
         required=True,
         type=str,
         action='append')

    parser.add_argument \
        ("--nprim",
         help="Number of primitive dihedrals",
         required=True,
         type=int)


    parser.add_argument \
        ("--opy",
         help="Output python script",
         required=True,
         type=str)


    
    parser.add_argument \
        ("--stride",
         help="Stride when reading scans from file. Default: 1",
         default=1,
         type=int)

    parser.add_argument \
        ("--nlmaxiter",
         help="Maximum number of nonlinear optimization steps. Default: 200",
         default=200,
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
    
    parser.add_argument \
        ("--no-overwrite-guess",
         help="If present, then the LL scan geometries are NOT replaced every time a new trial set of parameters are tested.",
         action='store_true')

    parser.add_argument \
        ("--bytype",
         help="Group parameters by atom type rather than atom index quartets",
         action='store_true')

    

    
    AddStandardOptions(parser)

    args = parser.parse_args()
    stdargs = StandardArgs(args)

    if stdargs.model != "sander":
        raise Exception("Fits can only be performed with sander")

    scans = []
    for llscan,hlscan in zip(args.llscan,args.hlscan):
        llgeoms = ase.io.read(llscan,format="extxyz",index=":")[::args.stride]
        hlgeoms = ase.io.read(hlscan,format="extxyz",index=":")[::args.stride]
        for x in llgeoms:
            x.info["energy"] = x.get_potential_energy()
        for x in hlgeoms:
            x.info["energy"] = x.get_potential_energy()
        pairs = []
        for llgeom,hlgeom in zip(llgeoms,hlgeoms):
            pairs.append(ScanPair(stdargs,llgeom,hlgeom))
        scans.append(pairs)
        
    udscans = MakeUniqueDihedralScans(stdargs,scans)
        
    NonlinearSolve(stdargs,udscans)                      
    WriteParmedScript(args.opy,stdargs.mol,
                      [ udscan.dfcn for udscan in udscans ],
                      args.bytype)

    
