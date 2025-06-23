#!/usr/bin/env python3

def FindDefaultValue(key,parser):
    default_value_of_my_arg = None
    for action in parser._actions:
        if action.dest == key:
            default_value_of_my_arg = action.default
            break
    return default_value_of_my_arg

def GetNondefaultArgs(args,parser,skip=[]):
    s = []
    for name,value in vars(args).items():
        dvalue = FindDefaultValue(name,parser)
        myname = name.replace("_","-")
        if myname in skip:
            continue
        if dvalue != value:
            if value != True:
                s.append( [f"--{myname}",f"{value}"] )
            else:
                s.append( [f"--{myname}"] )

    #return s
    return " ".join( [ " ".join(f) for f in s ] )




class Parameter(object):
    def __init__(self,mol,idxs):
        self.idxs=idxs
        self.res=mol.atoms[idxs[0]].residue.name
        self.names=[mol.atoms[i].name for i in idxs]
        self.types=[mol.atoms[i].type for i in idxs]
        self.instances = []
        
    def GetIdxStr(self):
        return "-".join(["%i"%(x) for x in self.idxs])

    def GetNameStr(self):
        return "-".join(self.names)

    def GetTypeStr(self):
        return "-".join(self.types)

    
    def GetParamByName(self):
        return "%s_%s"%(self.res,self.GetNameStr())
     
    def GetParamByType(self):
        return "%s_%s"%(self.res,self.GetTypeStr())


    
    def GetNameMasks(self):
        return [ f"@{name}" for name in self.names ]

    def GetTypeMasks(self):
        return [ f"@%%{name}" for name in self.types ]

    def AsDict(self):
        return { self.GetParamByType():
                 { "nprim": args.nprim,
                   "masks": None } }
    


if __name__ == "__main__":
    import sys
    import json
    import copy
    import argparse
    from ffpopt.Options import AddStandardOptions,StandardArgs
    
    parser = argparse.ArgumentParser \
        ( formatter_class=argparse.RawDescriptionHelpFormatter,
          description="""Perform sets up a workflow script for dihedral scans and fits""" )

    parser.add_argument\
        ("--bond",
         help="Two 0-based integers separated by a comma. This option can be used more than once",
         action='append',
         required=True)

    parser.add_argument\
        ("--delta",
         help="Scan step size. Default: 10 degrees",
         type=int,
         default=10)
    
    parser.add_argument\
        ("--nprim",
         help="Number of dihedral primitives. Default: 3",
         type=int,
         default=3)
    
    parser.add_argument\
        ("--maxiter",
         help="Number of training iterations. Each iteration repeats the sander scans. Default: 2",
         type=int,
         default=2)
    
    parser.add_argument\
        ("--bytype",
         help="If present, then",
         action='store_true')
    
    
    AddStandardOptions(parser)

    args = parser.parse_args()
    
    stdargs = StandardArgs(args)

    if args.model == "sander":
        raise Exception("Invalid argument: --model=sander\nChoose a different high-level model")
    
    useropts = GetNondefaultArgs(args,parser,
                                 skip=["bond","delta","bytype",
                                       "nprim","maxiter",
                                       "parm","crd","model"])

    bonds = [ [int(x) for x in bond.split(",") ] for bond in args.bond ]

    output = "global.frcmod"
    params = {}
    systems = []

    s = { "parm": None,
          "crd": args.crd,
          "output": None,
          "params": {},
          "profiles": [] }

    scans = []

    allparams = []
    ps = {}
    for bond in bonds:
        made_scan = False
        for d in stdargs.mol.dihedrals:
            if d.improper:
                continue
            idxs = [d.atom1.idx,d.atom2.idx,d.atom3.idx,d.atom4.idx]
            myidxs = []
            if idxs[1] == bond[0] and idxs[2] == bond[1]:
                myidxs = idxs
            elif idxs[2] == bond[0] and idxs[1] == bond[1]:
                myidxs = idxs[::-1]
            else:
                continue
            idxs = myidxs
            p = Parameter(stdargs.mol,idxs)
            name = p.GetParamByType()
            if name not in ps:
                ps[name] = p
            ps[name].instances.append( p.GetNameMasks() )
            allparams.append( name )

            if not made_scan:
                scans.append(p)
                made_scan=True

    uparams = list(set(allparams))
    if not args.bytype:
        for name in ps:
            s["params"][name] = ps[name].instances

        for name in uparams:
            params[name] = { "nprim": args.nprim, "masks": None }
    else:
        
        for name in uparams:
            typestr = name.split("_")[1]
            ts = [ f"@%{t}" for t in typestr.split("-") ]
            params[name] = { "nprim": args.nprim,
                             "masks": [ts]  }

    fh = sys.stdout

    fh.write("#!/bin/bash\n")
    fh.write("set -e\n")
    fh.write("set -u\n")

    for scan in scans:
        hlname = args.model.replace("/","_")
        oscan = hlname + "_" + scan.GetIdxStr()
        idxs = ",".join(["%i"%(x) for x in scan.idxs])
        opts = useropts

        fh.write("\n")
        fh.write(f"if [ ! -e {oscan}.xyz ]; then\n")
        fh.write(f"    echo Creating {oscan}.xyz\n")
        fh.write(f"    ffpopt-DihedScan.py -p {args.parm} -c {args.crd} --model {args.model} --dihed=\"{idxs}\" \\\n")
        fh.write(f"        --oscan {oscan}.xyz --delta {args.delta} {opts}\n\n")
        fh.write(f"fi\n")
        fh.write("\n")


    for scan in scans:
        oscan = "orig" +"_" + scan.GetIdxStr()
        idxs = ",".join(["%i"%(x) for x in scan.idxs])
        #opts = " ".join([ " ".join(useropts)])
        opts=useropts
        
        fh.write("\n")
        fh.write(f"if [ ! -e {oscan}.xyz ]; then\n")
        fh.write(f"    echo Creating {oscan}.xyz\n")
        fh.write(f"    ffpopt-DihedScan.py -p {args.parm} -c {args.crd} --model=sander --dihed=\"{idxs}\" \\\n")
        fh.write(f"        --oscan {oscan}.xyz --delta {args.delta} {opts}\n\n")
        fh.write(f"fi\n")
        fh.write("\n")
        
        
    
    
    for it in range(args.maxiter):
    
        
        pitname = "it%02i"%(it)
        citname = "it%02i"%(it+1)
        ss = copy.deepcopy(s)
        if it == 0:
            parm = args.parm
            #scanname = "orig"
        else:
            parm = pitname + ".parm7"
            #scanname = citname
        scanname = citname
        ss["output"] = citname + ".py"
        ss["parm"] = parm


        hlname = args.model.replace("/","_")

        ss["profiles"] = []
        for scan in scans:
            if it == 0:
                llscan = "orig"
            else:
                llscan = pitname
            prof = { "hl": "%s_%s.xyz"%(hlname,scan.GetIdxStr()),
                     "ll": "%s_%s.xyz"%(llscan,scan.GetIdxStr()),
                     "name": citname,
                     "plots": [ f"{scan.GetParamByType()}" ] }
            ss["profiles"].append(prof)
        
        datadict = { "params": params,
                     "output": f"{citname}.frcmod",
                     "systems": [ss] }

        jfh = open(f"{citname}.json","w")
        json.dump(datadict,jfh,indent=4)
        jfh.close()

        fh.write("\n")
        fh.write(f"if [ ! -e {citname}.py ]; then\n")
        fh.write(f"    echo Creating {citname}.py\n")
        fh.write(f"    ffpopt-GenDihedFit.py --ase-opt {citname}.json\n")
        fh.write("fi\n")
        fh.write("\n")

        fh.write("\n")
        fh.write(f"if [ ! -e {citname}.parm7 ]; then\n")
        fh.write(f"    echo Creating {citname}.parm7\n")
        fh.write(f"    python3 {citname}.py {args.parm} {citname}.parm7\n")
        fh.write("fi\n")
        fh.write("\n")

        
        for scan in scans:
            
            oscan = scanname +"_" + scan.GetIdxStr()
            idxs = ",".join(["%i"%(x) for x in scan.idxs])
            #opts = " ".join([ " ".join(useropts)])
            opts=useropts
            fh.write("\n")
            fh.write(f"if [ ! -e {oscan}.xyz ]; then\n")
            fh.write(f"    echo Creating {oscan}.xyz\n")
            fh.write(f"    ffpopt-DihedScan.py -p {citname}.parm7 -c {args.crd} --model=sander --dihed=\"{idxs}\" \\\n")
            fh.write(f"        --oscan {oscan}.xyz --delta {args.delta} {opts}\n\n")
            fh.write(f"fi\n")
            fh.write("\n")

    fh.close()
