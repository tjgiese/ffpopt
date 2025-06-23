#!/usr/bin/env python3

class Constraint(object):

    def __init__(self,idxs,value=None,graph=None):
        self.idxs = idxs
        self.value = value
        self.mask = None
        if graph is not None:
            self.SetMask(graph)

    @classmethod
    def from_str(cls,s,graph=None):
        value = None
        if "=" in s:
            s,value = s.split("=")
            value = float(value)
        idxs = [int(x) for x in s.split(",")]
        return cls(idxs,value,graph=graph)

    def SetMask(self,graph):
        from . AmberParm import RotateMask
        self.mask = None
        if len(self.idxs) == 4:
            self.mask = RotateMask(graph,self.idxs)
            
    def __str__(self):
        s = ",".join(["%i"%(x) for x in self.idxs])
        if self.value is not None:
            s += f"={self.value}"
        return s

    def fill(self,atoms,force=False):
        import copy
        out = copy.deepcopy(self)
        if self.value is None or force:
            if len(self.idxs) == 2:
                out.value = atoms.get_distance(self.idxs[0],self.idxs[1])
            elif len(self.idxs) == 3:
                out.value = atoms.get_angle(self.idxs[0],self.idxs[1],self.idxs[2])
            elif len(self.idxs) == 4:
                out.value = atoms.get_dihedral(self.idxs[0],self.idxs[1],self.idxs[2],self.idxs[3])
        return out

    def modify(self,atoms):
        out = atoms.copy()
        if self.value is not None:
            if len(self.idxs) == 2:
                out.set_distance(self.idxs[0],self.idxs[1],self.value)
            elif len(self.idxs) == 3:
                out.set_angle(self.idxs[0],self.idxs[1],self.idxs[2],self.value)
            elif len(self.idxs) == 4:
                out.set_dihedral(self.idxs[0],self.idxs[1],self.idxs[2],self.idxs[3],self.value,mask=self.mask)
        return out

    def to_geometric(self):
        if self.value is None:
            raise Exception("self.value should not be None. Use the fill method to assign values")
        s = " ".join( [f"{1+x}" for x in self.idxs] )
        if len(self.idxs) == 2:
            s = "distance "+s+f" {self.value}"
        elif len(self.idxs) == 3:
            s = "angle "+s+f" {self.value}"
        elif len(self.idxs) == 4:
            s = "dihedral "+s+f" {self.value}"
        return s


def ApplyConstraints(atoms,cons,graph=None):
    out = atoms.copy()
    for con in cons:
        if graph is not None:
            out.SetMask(graph)
        out = con.modify(out)
    return out


def FillConstraints(atoms,cons,force=False):
    return [ con.fill(atoms,force=force) for con in cons ]
        

def constraints2geometric(cons):
    return [ con.to_geometric() for con in cons ]


def constraints2ase(cons):
    from ase.constraints import FixInternals
    bonds=[]
    angles=[]
    diheds=[]
    for con in cons:
        if len(con.idxs) == 2:
            bonds.append( [con.value,con.idxs] )
        elif len(con.idxs) == 3:
            angles.append( [con.value,con.idxs] )
        elif len(con.idxs) == 4:
            diheds.append( [con.value,con.idxs] )
    return FixInternals(bonds=bonds,angles_deg=angles,dihedrals_deg=diheds)


def constraints2info(atoms,cons):
    if cons is not None:
        atoms.info["constraints"] =  "[" + ",".join( ["\"%s\""%(con) for con in cons] ) + "]"


def info2constraints(atoms,graph=None):
    import json
    cons = None
    if "constraints" in atoms.info:
        cs = json.loads( atoms.info["constraints"] )
        cons = [ Constraint.from_str( c, graph=graph ) for c in cs ]
    return cons

