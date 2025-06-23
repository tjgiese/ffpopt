#!/usr/bin/env python3



def DeleteDihedrals(p,list_of_idxs):
    from parmed.tools.actions import deleteDihedral, addDihedral
    for idxs in list_of_idxs:
        cmd = deleteDihedral(p,
                             f"@{idxs[0]+1}",
                             f"@{idxs[1]+1}",
                             f"@{idxs[2]+1}",
                             f"@{idxs[3]+1}")
    
        cmd.execute()
        
    
    

def ChangeDihedrals(p,idxs,xs,fc=None,bytype=False):
    from collections import defaultdict as ddict
    from parmed.tools.actions import deleteDihedral, addDihedral

    if bytype:

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
        allseles = [ [ f"@{idx[i]+1}" for i in range(4) ]
                     for idx in allidxs ]

    else:
        allseles = [ [ f"@{i+1}" for i in idxs ] ]

        
    for seles in allseles:
        
        cmd = deleteDihedral(p,seles[0],seles[1],seles[2],seles[3])    
        cmd.execute()

        for x in xs:
            per = x.per
            phase = x.phase
            scee = x.scee
            scnb = x.scnb
            
            k = x.phi_k
            if fc is not None:
                k = fc
    
            cmd = addDihedral\
                (p,
                 seles[0],seles[1],seles[2],seles[3],
                 k,per,phase,scee,scnb)
            cmd.execute()

    

    

class PrimDihedFcn(object):
    def __init__(self,fc,phase,per):
        self.fc = fc
        self.phase = phase
        self.per = per

    def __str__(self):
        return f"({self.fc},{self.phase},{self.per})"
        
    def CptEne(self,ang):
        return self.fc * self.CptEterm(ang)

    def CptEterm(self,ang):
        import numpy as np
        a = (self.per * ang + self.phase)*(np.pi/180)
        return 1+np.cos(a)

    
class MultiDihedFcn(object):
    def __init__(self,idxs,prims):
        import copy
        self.idxs = copy.deepcopy(idxs)
        self.prims = copy.deepcopy(prims)

    def CptEne(self,ang):
        import numpy as np
        outs = []
        for p in self.prims:
            outs.append( p.CptEne(ang) )
        return np.sum(outs,axis=0)

    def SetFCs(self,fcs):
        if len(fcs) != len(self.prims):
            raise Exception(f"Size mismatch {len(fcs)} {len(self.prims)}")
        for i in range(len(fcs)):
            self.prims[i].fc = fcs[i]

    def ResetFCs(self):
        self.SetFCs( [1]*len(self.prims) )

    def __str__(self):
        return "[%s]"%(",".join([str(x) for x in self.prims]))


    
def CptDihedralEne(atoms,mdfs):
    import numpy as np
    enes = []
    for mdf in mdfs:
        idxs = mdf.idxs
        ang = atoms.get_dihedral(idxs[0],idxs[1],idxs[2],idxs[3])
        e = mdf.CptEne(ang)
        enes.append(e)
    return np.array(enes)
        



def FindDihedrals(p,idxs,impropers=False):
    tkey = tuple(idxs)
    xs=[]
    for x in p.dihedrals:
        if x.improper and not impropers:
            continue
        fkey = (x.atom1.idx,x.atom2.idx,x.atom3.idx,x.atom4.idx)
        rkey = tuple(list(fkey)[::-1])
        #print(tkey,fkey,rkey,fkey == tkey, rkey == tkey, x.atom1.type, x.atom2.type,x.atom3.type,x.atom4.type)
        if fkey == tkey or rkey == tkey:
            xs.append(x)
    return xs




def GetMultiDihedFcnFromIdxs(p,idxs):
    xs = FindDihedrals(p,idxs)
    prims = []
    for x in xs:
        t = x.type
        prims.append( PrimDihedFcn(t.phi_k,t.phase,t.per) )
    return MultiDihedFcn(idxs,prims)


def ChangeParmFromMultiDihedFcn(p,fcn):
    out = CopyParm(p)
    scee = 1.2
    scnb = 2.0
    xs = []
    for prim in fcn.prim:
        per = prim.per
        ph = prim.phase
        fc = prim.fc
        x = DihedralType(fc, per, ph, scee, scnb)
        xs.append(x)
    ChangeDihedrals(out,fcn.idxs,xs)
    return out




# def GetDihedClasses(idxs=None):
    
#     class1 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,3)] ) ]
               
#     class2 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,  0,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,180,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,180,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1, 0,2)] ) ]

#     class3 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,  0,2),
#                           PrimDihedFcn(1,  0,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,180,2),
#                           PrimDihedFcn(1,180,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,180,2),
#                           PrimDihedFcn(1,180,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,  0,2),
#                           PrimDihedFcn(1,180,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,180,2),
#                           PrimDihedFcn(1,  0,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,  0,2),
#                           PrimDihedFcn(1,180,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,180,2),
#                           PrimDihedFcn(1,  0,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,  0,2),
#                           PrimDihedFcn(1,  0,3)] ) ]


#     return { 1: class1,
#              2: class2,
#              3: class3 }




def GetDihedClasses(idxs=None):
    
    class1 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1)] ),
              MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,2)] ),
              MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,3)] ) ]
               
    # class2 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
    #                                 PrimDihedFcn(1,  0,2)] ) ]

    # class3 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
    #                                 PrimDihedFcn(1,  0,2),
    #                                 PrimDihedFcn(1,  0,3)] ) ]

    cs = {}
    cs[1] = class1
    for j in range(2,7):
        cs[j] = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,n) for n in range(1,j+1)])]
    

    return cs





def EnergyScansWithoutDihedrals(stdargs,list_of_scans,cons):
    from . AmberParm import CopyParm
    from . Dihedrals import DeleteDihedrals
    from . Dihedrals import GetMultiDihedFcnFromIdxs
    
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



def IsolatedLinearSolve(stdargs,idxs,llgeoms,hlenes,nprim,pname):
    #from . AmberParm import GetDihedClasses
    from . Constraints import FillConstraints
    from . Constraints import Constraint
    import numpy as np
    import copy
    from . constants import AU_PER_KCAL_PER_MOL
    from . constants import AU_PER_ELECTRON_VOLT

    KCAL_PER_EV = AU_PER_ELECTRON_VOLT() / AU_PER_KCAL_PER_MOL()
    
    dfcns = GetDihedClasses(idxs=idxs)[nprim]
    list_of_scans = [ llgeoms ]
    cons = [ Constraint(idxs,graph=stdargs.graph) ]
    list_of_enes = EnergyScansWithoutDihedrals(stdargs,list_of_scans,cons)
    #print(list_of_enes[0])
    llenes = np.array(list_of_enes[0])
    hlenes = np.array(hlenes,copy=True)

    llenes *= KCAL_PER_EV
    hlenes *= KCAL_PER_EV

    hlenes -= np.amin(hlenes)
    llenes -= np.amin(llenes)
    

    angs = []
    for igeom,llgeom in enumerate(llgeoms):
        o = FillConstraints(llgeom,cons,force=True)
        v = o[0].value
        angs.append(v)

    data = []
    for i in range(len(llgeoms)):
        data.append( [angs[i],llenes[i],hlenes[i]] )
    data = sorted(data,key=lambda x: x[0])
    angs = np.array( [x[0] for x in data] )
    llenes = np.array( [x[1] for x in data] )
    hlenes = np.array( [x[2] for x in data] )
    y = hlenes-llenes
    
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

    fh = open(f"iso.{pname}.dat","w")
    fh.write("# %s\n"%(str(bestdfcn)))
    for i in range(npts):
        fh.write("%12.3f %20.10e %20.10e %20.10e\n"%\
                 ( angs[i], hlenes[i], llenes[i],
                   llenes[i]+bestvalues[i] ) )
    fh.close()
    
    return bestdfcn
    


def AngularStdDev(angs):
    import numpy as np
    from scipy.stats import circstd
    rads = np.deg2rad(np.array(angs))
    return np.rad2deg(circstd(rads))



class ParamType(object):
    
    def __init__(self,name,nprim,masks):
        self.name = name
        self.nprim = nprim
        self.masks = masks
        self.dfcns = None


        
    def get_dict(self):
        d = {}
        d[self.name] = { "nprim": self.nprim,
                         "masks": self.masks }
        return d

    
        
class ParamInstance(object):
    
    def __init__(self,ptype,masks,mol,impropers=False):
        from parmed.amber.mask import AmberMask
        
        self.ptype = ptype
        self.masks = masks
        self.dihedidxs = []
        
        possible_diheds = []
        for mask in masks:
            ats = []
            for i in range(4):
                sidxs = [idx for idx in AmberMask(mol,mask[i]).Selected()]
                ats.append(sidxs)
            for aidx in ats[0]:
                #amask = f"@{aidx+1}"
                for bidx in ats[1]:
                    #bmask = f"@{bidx+1}"
                    for cidx in ats[2]:
                        #cmask = f"@{cidx+1}"
                        for didx in ats[3]:
                            #dmask = f"@{didx+1}"
                            fwd = (aidx,bidx,cidx,didx)
                            possible_diheds.append(fwd)
            for tgt in possible_diheds:
                for d in mol.dihedrals:
                    if d.improper and not impropers:
                        continue
                    fwd = (d.atom1.idx,d.atom2.idx,d.atom3.idx,d.atom4.idx)
                    rev = tuple(list(fwd)[::-1])
                    if fwd == tgt or rev == tgt:
                        if tgt not in self.dihedidxs:
                            self.dihedidxs.append(tgt)
                            


    def get_dict(self):
        from collections import defaultdict as ddict
        d = ddict(list)
        name = self.ptype.name
        for didx in self.dihedidxs:
            d[name].append( [ f"@{didx[0]+1}",f"@{didx[1]+1}",
                              f"@{didx[2]+1}",f"@{didx[3]+1}" ] )
        return d


    
class ProfileType(object):
    def __init__(self,hl,ll,name,plots,stdargs):
        from . Reader import ReadGeomsFromXYZ
        self.hl = hl
        self.ll = ll
        self.name = name
        self.plots = plots
        self.hlatoms,self.hlcons = ReadGeomsFromXYZ(stdargs,self.hl)
        self.llatoms,self.llcons = ReadGeomsFromXYZ(stdargs,self.ll)
        
        if len(self.hlatoms) != len(self.llatoms):
            raise Exception(f"Structure count mismatch in {self.hl} "
                            +f"and {self.ll} ({len(self.hlatoms)} vs "
                            +f"{len(self.llatoms)})")
        
        s = stdargs.args.stride
        self.hlatoms = self.hlatoms[::s]
        self.hlcons = self.hlcons[::s]
        self.llatoms = self.llatoms[::s]
        self.llcons = self.llcons[::s]


        
    def get_dict(self):
        return { "hl": self.hl, "ll": self.ll,
                 "name": self.name, "plots": self.plots }

        
        
class SystemType(object):
    
    def __init__(self,stdargs,output,profilelist,pdict,ptypedict):
        self.stdargs = stdargs
        self.output = output
        self.pinstances = []
        self.profiles = []
        
        for pname in pdict:
            if pname not in ptypedict:    
                raise Exception(f"Parameter {pname} in {stdargs.parm} "
                                +f"not found in parameter list "
                                +f"{[name for name in ptypedict]}")
            
            ptype = ptypedict[pname]
            masks = pdict[pname]
            self.pinstances.append( ParamInstance(ptype,masks,stdargs.mol) )

        for pname in ptypedict:
            if pname not in pdict:
                ptype = ptypedict[pname]
                masks = ptype.masks
                pinst = ParamInstance(ptype,masks,stdargs.mol)
                if len(pinst.dihedidxs) > 0:
                    self.pinstances.append(pinst)
            
        for profile in profilelist:
            for key in ["hl","ll","name","plots"]:
                if key not in profile:
                    raise Exception(f"'profiles' section is missing '{key}' key")
            hl = profile["hl"]
            ll = profile["ll"]
            name = profile["name"]
            plots = profile["plots"]
            self.profiles.append( ProfileType(hl,ll,name,plots,self.stdargs) )

            

    def make_new_parm(self):
        from parmed import DihedralType
        from . Dihedrals import ChangeDihedrals
        from . AmberParm import CopyParm

        p = CopyParm(self.stdargs.mol)
        
        scee = 1.2
        scnb = 2.0
        for pinst in self.pinstances:
            xs = []
            for prim in pinst.ptype.dfcns.prims:
                per = prim.per
                ph = prim.phase
                fc = prim.fc
                x = DihedralType(fc, per, ph, scee, scnb)
                xs.append(x)
            for idxs in pinst.dihedidxs:
                ChangeDihedrals(p,idxs,xs)
        return p

    
            
    def get_dict(self):
        params = {}
        profiles = []

        for p in self.pinstances:
            x = p.get_dict()
            for key in x:
                params[key] = x[key]

        for p in self.profiles:
            profiles.append(p.get_dict())
        
        d = { "parm": self.stdargs.parm,
              "crd": self.stdargs.crd,
              "output": self.output,
              "params": params,
              "profiles": profiles }
        return d


    
    def find_pinstance(self,pname):
        p=None
        for pinst in self.pinstances:
            if pinst.ptype.name == pname:
                p = pinst
        return p


    def write_output(self):
        import copy
        dfcns = []
        for pinst in self.pinstances:
            for idxs in pinst.dihedidxs:
                dfcn = copy.deepcopy(pinst.ptype.dfcns)
                dfcn.idxs = idxs
                dfcns.append(dfcn)
        WriteParmedScript(self.output,self.stdargs.mol,dfcns)
    
    
class FitInputType(object):
    
    @classmethod
    def from_file(cls,args,fname):
        import json
        data = {}
        with open(fname,"r") as file:
            data = json.load(file)
        return cls(args,data)
    
        
    def __init__(self,args,datadict):
        from . Options import StandardArgs
        self.iteration = 0

        for key in ["params","output","systems"]:
            if key not in datadict:
                raise Exception(f"Input is missing '{key}' key")
                
        self.ptypedict = {}
        for pname in datadict["params"]:
            #print(datadict["params"][pname])
            nprim = datadict["params"][pname]["nprim"]
            masks = datadict["params"][pname]["masks"]
            if masks is not None:
                for mask in masks:
                    for atmask in mask:
                        if atmask[:2] != "@%":
                            raise Exception(f"Global dihedral {mask} contains mask {atmask}, which should start with @%%")
            self.ptypedict[pname] = ParamType(pname,nprim,masks)
        self.output = datadict["output"]
        self.systems = []
        for s in datadict["systems"]:
            for key in ["parm","crd","output","params","profiles"]:
                if key not in s:
                    raise Exception(f"'systems' section is missing '{key}' key")
            parm = s["parm"]
            crd = s["crd"]
            output = s["output"]
            pdict = s["params"]
            profilelist = s["profiles"]
            stdargs = StandardArgs.from_manual_parm(parm,crd,args)
            self.systems.append\
                ( SystemType\
                  (stdargs,output,profilelist,pdict,self.ptypedict) )

        unused = []
        for pname in self.ptypedict:
            found=False
            for s in self.systems:
                for pinst in s.pinstances:
                    if pinst.ptype.name == pname:
                        found=True
            if not found:
                unused.append(pname)
        for pname in unused:
            print(f"Removing parameter {pname} because it was unused")
            del self.ptypedict[pname]
                
            
    def get_dict(self):
        d = {}
        d["params"] = {}
        for pname in self.ptypedict:
            x = self.ptypedict[pname].get_dict()
            d["params"][pname] = x[pname]
        d["output"] = self.output
        d["systems"] = []
        for s in self.systems:
            d["systems"].append( s.get_dict() )

        return d

    def get_json(self):
        import json
        d = self.get_dict()
        return json.dumps(d,indent=4)
            

    def make_initial_guesses(self):
        for pname in self.ptypedict:
            #print(pname)
            bests = None
            bestpinst = None
            bestprof = None
            bestidxs = None
            beststd = -1
            
            for s in self.systems:
                for pinst in s.pinstances:
                    if pinst.ptype.name == pname:
                        for idxs in pinst.dihedidxs:
                            for prof in s.profiles:
                                angs = []
                                for atoms in prof.llatoms:
                                    ang = atoms.get_dihedral(*idxs)
                                    angs.append(ang)
                                    #print(angs)
                                if len(angs) > 2:
                                    astd = AngularStdDev(angs)
                                    #print(astd,beststd)
                                    if astd > beststd:
                                        beststd = astd
                                        bestprof = prof
                                        bestpinst = pinst
                                        bestidxs = idxs
                                        bests = s
            llgeoms = bestprof.llatoms
            hlenes = [ hlgeom.info["energy"]
                       for hlgeom in bestprof.hlatoms ]
            nprim = bestpinst.ptype.nprim
            
            dfcns = IsolatedLinearSolve\
                (bests.stdargs,bestidxs,llgeoms,hlenes,nprim,pname)
            
            bestpinst.ptype.dfcns = dfcns
        return self.get_params()

            
    def get_num_params(self):
        n = 0
        for pname in self.ptypedict:
            n += self.ptypedict[pname].nprim
        return n

    

    def get_params(self):
        import numpy as np
        x=[]
        for pname in self.ptypedict:
            ptype = self.ptypedict[pname]
            dfcns = ptype.dfcns
            for prim in dfcns.prims:
                x.append( prim.fc )
        return np.array(x)


    
    def set_params(self,x):
        ipar = 0
        for pname in self.ptypedict:
            ptype = self.ptypedict[pname]
            nprim = ptype.nprim
            ptype.dfcns.SetFCs( x[ipar:ipar+nprim] )
            ipar += nprim

    def write_output(self):
        from parmed.amber import AmberParameterSet
        from parmed import DihedralTypeList
        from parmed import DihedralType
        #from parmed.amber.mask import AmberMask
        
        for s in self.systems:
            s.write_output()
            
        fmod = AmberParameterSet()
        for pname in self.ptypedict:
            ptype = self.ptypedict[pname]
            if ptype.masks is not None and ptype.dfcns is not None:

                typs = DihedralTypeList()
                scee = 1.2
                scnb = 2.0
                for iprim,prim in enumerate(ptype.dfcns.prims):
                    per = prim.per
                    ph = prim.phase
                    fc = prim.fc
                    typ = DihedralType(fc, per, ph, scee, scnb)
                    typs.append(typ)
                
                for mask in ptype.masks:
                    atypes = [ x[2:] for x in mask ]
                    fwd = tuple(atypes)
                    rev = tuple(list(fwd)[::-1])
                    fmod.dihedral_types[fwd] = typs
                    fmod.dihedral_types[rev] = typs
                    
        fmod.write(self.output)
    

def DihedFitObjFcn(x,self):
    import numpy as np
    from . GeomOpt import GeomOpt
    from . constants import AU_PER_KCAL_PER_MOL
    from . constants import AU_PER_ELECTRON_VOLT

    KCAL_PER_EV = AU_PER_ELECTRON_VOLT() / AU_PER_KCAL_PER_MOL()
  
    
    chisq = 0

    it = self.iteration
    
    self.set_params(x)
    for isys,s in enumerate(self.systems):
        p = s.make_new_parm()
        p.save("tmp.parm7",overwrite=True)
        s.stdargs.calc = s.stdargs.MakeCalc(parm="tmp.parm7")

        for iprof,prof in enumerate(s.profiles):
            llene = []
            hlene = []
            for igeom in range(len(prof.llatoms)):
                atoms = prof.llatoms[igeom].copy()
                cons = prof.llcons[igeom]
                atoms = GeomOpt(atoms,s.stdargs,constraints=cons)
                hlene.append( prof.hlatoms[igeom].info["energy"] * KCAL_PER_EV )
                llene.append( atoms.info["energy"] * KCAL_PER_EV )
                prof.llatoms[igeom].set_positions( atoms.get_positions() )
                
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

            if prof.name is None or prof.plots is None:
                continue
            elif len(prof.plots) == 0:
                continue

            for pname in prof.plots:
                pinst = s.find_pinstance(pname)
                if pinst is None:
                    continue
                            
                for idxs in pinst.dihedidxs:
                    angs = []
                    for atoms in prof.llatoms:
                        ang = atoms.get_dihedral( *idxs )
                        angs.append(ang)
                    data = []
                    for i in range(len(angs)):
                        data.append( [angs[i],hlene[i],llene[i]] )
                    data = sorted(data,key=lambda x: x[0])
                
                    idxsname = "-".join( [ f"{i}" for i in idxs] )
                    fname = f"mfit.{prof.name}.{idxsname}.{it:04d}.dat"
                    fh = open(fname,"w")
                    fh.write("# %25.14f\n"%(mychisq))
                    for x in data:
                        fh.write("%20.10e %20.10e %20.10e\n"%(x[0],x[1],x[2]))
                    fh.close()
                
    self.iteration += 1
    return chisq

                


def NonlinearSolve(args,finp):
    from scipy.optimize import minimize
    import numpy as np
    import copy

    #objfcn = NonlinearObjective(stdargs,udscans)
    
    n = finp.get_num_params()
    x = finp.make_initial_guesses()
    xlo = x[:] - 2.
    xhi = x[:] + 5.
    bounds = [ (lo,hi) for lo,hi in zip(xlo,xhi) ]
    
    res = minimize( DihedFitObjFcn, x, args=(finp,),
                    method='COBYLA',
                    bounds=bounds,
                    options={ "rhobeg": args.nlrhobeg,
                              "tol": args.nltol,
                              "maxiter": args.nlmaxiter,
                              "disp": True })

    print(res)

    finp.set_params(res.x)



def WriteParmedScript(fname,p,dfcns): #,bytype):
    from collections import defaultdict as ddict
    from . Dihedrals import FindDihedrals
    
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
    

    #if not bytype:
    if True:
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

    
