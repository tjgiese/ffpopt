#!/usr/bin/env python3

from collections import defaultdict as ddict

def parmed2ase(mol):
    from parmed import periodic_table
    import numpy as np
    import ase
    charge = int(round(sum([ a.charge for a in mol.atoms ])))
    eles = [ periodic_table.Element[a.element]
            for a in mol.atoms]
    crds = np.array([ [ a.xx, a.xy, a.xz ] for a in mol.atoms ])
    atlist = "".join( ["%s1"%(ele) for ele in eles ] )
    return ase.Atoms(atlist,positions=crds),charge
    


def parmed2graph(mol):
    import parmed
    import numpy as np
    edges = []
    for x in mol.bonds:
        edges.append( "%i~%i"%(x.atom1.idx,x.atom2.idx) )
    return GraphSearch(edges)
    

def CopyParm( parm ):
    import copy
    try:
        parm.remake_parm()
    except:
        pass
    p = copy.copy( parm )
    p.coordinates = copy.copy( parm.coordinates )
    p.box = copy.copy( parm.box )
    try:
        p.hasbox = copy.copy( parm.hasbox )
    except:
        p.hasbox = False
    return p


def RotateMask(graph,idxs):
    left = "%i"%(idxs[1])
    right = "%i"%(idxs[2])
    nat = len(graph.nodes)
    gleft = []
    gright = []
    for n in graph.nodes:
        if n == left:
            gleft.append(int(n))
        elif n == right:
            gright.append(int(n))
        else:
            rleft = len(graph.FindMinPaths(left,n)[0])
            rright = len(graph.FindMinPaths(right,n)[0])
            if rleft < rright:
                gleft.append(int(n))
            else:
                gright.append(int(n))
    mask = [0]*nat
    if len(gleft) < len(gright):
        gmove = gleft
    else:
        gmove = gright
    for i in gmove:
        mask[i] = 1

    if mask[idxs[3]] == 0:
        mask = [ 1-x for x in mask ]
        
    return mask



class GraphSearch(object):
    """A class that performs path traversals of a graph

    Attributes
    ----------
    nodes : list of str
        A name for each unique node in the graph

    edges : dict of list
        The graph edges. The keys are the edge name; e.g., "a~b",
        and the values are the list of nodes in the edge; e.g., ["a","b"]

    Methods
    -------
    """
    
    def __init__(self,edges):
        from collections import defaultdict as ddict
        self.nodes = list(set([ x for e in edges for x in e.split("~") ]))
        self.nodes.sort()
        self.edges = ddict(list)
        for e in edges:
            a,b=e.split("~")
            self.edges[a].append(b)
            self.edges[b].append(a)

            
    def PathToStr(self,path):
        """Return a path name, given a path

        Parameters
        ----------
        path : list of str
            The list of node names defining a path; e.g., ["a","b"]

        Returns
        name : str
            The path name; e.g., "a~b"
        """
        return "~".join(path)

    
    def StrToPath(self,name):
        """Return a path, given its name

        Parameters
        ----------
        name : str
            The path name; e.g., "a~b"

        Returns
        -------
        path : list of str
            The list of node names defining a path; e.g., ["a","b"]
        """
        return name.split("~")

    
    def FindAllPaths(self,snode,tnode,minsize=2):
        """Return a list of all paths that connect snode to tnode

        Parameters
        ----------
        snode : str
            The name of the starting node

        tnode : str
            The name of the ending node

        minsize : int, default=2
            The paths that reference fewer than minsize nodes are
            excluded from the returned list of paths

        Returns
        -------
        paths : list of list of str
            The found paths. Each path is a list of nodes
        """
        from collections import defaultdict as ddict
        if snode not in self.nodes:
            raise Exception("Node %s is not in graph"%(snode))
        if tnode not in self.nodes:
            raise Exception("Node %s is not in graph"%(tnode))
        allpaths = self._find(snode,tnode,path=[],
                              visited=ddict(lambda: False),
                              allpaths=[])
        if len(allpaths) > 0:
            if minsize > 2:
                sizes = np.array([ len(path) for path in allpaths ],dtype=int)
                allpaths = [allpaths[x] for x in np.where( sizes >= minsize )[0]]
            allpaths.sort()
        return allpaths

    
    def FindMinPaths(self,snode,tnode,minsize=2):
        """Return a list of all minimum-length paths that 
        connect snode to tnode. Although there are many long and short
        pathways that may connect the nodes, this function only returns
        the list of paths that reference the fewest possible number of
        nodes to connect the endpoints.

        Parameters
        ----------
        snode : str
            The name of the starting node

        tnode : str
            The name of the ending node

        minsize : int, default=2
            The returned minimum-length paths will not be smaller than
            minsize

        Returns
        -------
        paths : list of list of str
            The found paths. Each path is a list of nodes
        """
        import numpy as np
        allpaths = [ x for x in self.FindAllPaths(snode,tnode)
                     if len(x) >= minsize ]
        if len(allpaths) > 0:
            sizes = np.array([ len(path) for path in allpaths ],dtype=int)
            return [allpaths[x]
                    for x in np.where( sizes == sizes.min() )[0]]
        else:
            return []
    
    def FindAllCycles(self,minsize=3):
        """Returns all unique closed cycles within a graph

        Parameters
        ----------
        minsize : int, default=3
            All cycles that contain fewer than minsize nodes will be
            excluded

        Returns
        -------
        paths : list of list of str
            The found paths. Each path is a list of nodes
        """
        import numpy as np
        cycs=[]
        for snode in self.nodes:
            for tnode in self.edges[snode]:
                paths = self.FindAllPaths(snode,tnode,minsize=3)
                for path in paths:
                    minnode = min(path)
                    minidx = [i for i,j in enumerate(path) if j == minnode][0]
                    path = path[minidx:] + path[:minidx]
                    if path[-1] < path[1]:
                        f = path.pop(0)
                        path.append(f)
                        path.reverse()
                    path.append(path[0])
                    cycs.append( self.PathToStr(path) )
        allpaths = [ self.StrToPath(x) for x in list(set(cycs)) ]
        if len(allpaths) > 0:
            if minsize > 2:
                sizes = np.array([ len(path) for path in allpaths ],dtype=int)
                allpaths = [allpaths[x] for x in np.where( sizes >= minsize )[0]]
            allpaths.sort()            
        return allpaths

    
    def FindMinCycles(self):
        """Returns all minimum-length unique closed cycles within a graph.
        That is, the returned cycles cannot be expressed as a sum of two
        smaller cycles.

        Parameters
        ----------
        minsize : int, default=3
            All cycles that contain fewer than minsize nodes will be
            excluded

        Returns
        -------
        paths : list of list of str
            The found paths. Each path is a list of nodes
        """
        cycs=[]
        for snode in self.nodes:
            for tnode in self.edges[snode]:
                paths = self.FindMinPaths(snode,tnode,minsize=3)
                for path in paths:
                    minnode = min(path)
                    minidx = [i for i,j in enumerate(path) if j == minnode][0]
                    path = path[minidx:] + path[:minidx]
                    if path[-1] < path[1]:
                        f = path.pop(0)
                        path.append(f)
                        path.reverse()
                    path.append(path[0])
                    cycs.append( self.PathToStr(path) )
        allpaths = [ self.StrToPath(x) for x in list(set(cycs)) ]
        allpaths.sort()
        return allpaths

                
    def _find(self,start,stop,path=[],
              visited=ddict(lambda: False),
              allpaths=[]):
        """Utility function that recursively traverses the graph to
        find paths.  One should instead use the FindAllPaths method,
        which uses this recursive algorithm, but which is harder to
        screw up, because of python's weird behavior of saved-state
        in recursive functions

        Parameters
        ----------
        start : str
            The starting node

        stop : str
            The stopping node

        path : list of str
            The current list of nodes in the path 
            (should be set as [] when calling)

        visited : dict of bool
            Indicates if a node has already been used within the path

        allpaths : list of list of str
            The list of all paths found in the graph
            (should be set to [] when calling)
        
        Returns
        -------
        allpaths : list of list of str
            The list of all paths found in the graph
        """

        #
        # Adapted from
        # https://www.geeksforgeeks.org/find-paths-given-source-destination
        # Recursive Depth First Traversal with a boolean array used to avoid
        # revisiting a node twice
        #
        visited[start]=True
        path.append(start)
        if start == stop:
            q = [ p for p in path ]
            n = len(q)
            skip=False
            for h in allpaths:
                if len(h) == n:
                    same=True
                    for i in range(n):
                        if q[i] != h[i]:
                            same=False
                            break
                    if same:
                        skip=True
            if not skip:
                allpaths.append( q )
        else:
            for i in self.edges[start]:
                if not visited[i]:
                    allpaths=self._find(i,stop,path,visited,allpaths=allpaths)
        path.pop()
        visited[start]=False
        return allpaths

