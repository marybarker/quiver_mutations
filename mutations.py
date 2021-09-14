import numpy as np
import copy
import networkx as nx
import matplotlib.pyplot as plt

def flattenList(l):
    if type(l) is list or type(l) is tuple:
        return [x for y in l for x in flattenList(y)]
    return [l]


class QuiverWithPotential():

    def __init__(self, graph_obj, potential=None, positions=None):
        if isinstance(graph_obj, list): # if graph input is a list of edges
            self.Q1 = graph_obj
            self.incidence_matrix = matrixFromEdges(graph_obj)
        elif isinstance(graph_obj, np.matrix): # if graph input is an incidence matrix
            self.incidence_matrix = graph_obj
            self.Q1 = edgesFromMatrix(graph_obj)
        else:
            self.incidence_matrix = graph_obj.incidence_matrix
            self.Q1 = graph_obj.Q1
        self.positions=positions
        self.Q0 = range(self.incidence_matrix.shape[0])

        self.potential = {}
        if potential is not None:
            if isinstance(potential, dict):
                self.potential = {cycleOrder(tuple(k)):v for k,v in potential.items()}
            elif isinstance(potential, list) and (len(potential) == 2):
                self.add_term_to_potential(potential[0], potential[1])

        self.arrows_with_head = [[] for x in self.Q0]
        self.arrows_with_tail = [[] for x in self.Q0]
        self.loops_at = [[] for x in self.Q0]
        self.can_mutate = [True for x in self.Q0]

        for i, e in enumerate(self.Q1):
            if e[0] != e[1]:
                self.arrows_with_head[e[1]].append((i,e))
                self.arrows_with_tail[e[0]].append((i,e))
            else:
                self.loops_at[e[0]].append((i,e))
                self.can_mutate[e[0]] = False


    def print_potential(self):
        return " + ".join(["%d X_"%v+".".join([str(self.Q1[x][0]) for x in k]) for k,v in self.potential.items()])


    def add_term_to_potential(self, edge_cycle, coef=1, input_format="vertex_order"):
        if input_format == "vertex_order":
            if not isinstance(edge_cycle[0], tuple):
                edge_cycle = [edge_cycle]

            if coef == 1:
                coef = [1 for e in edge_cycle]
            if len(coef) != len(edge_cycle):
                print("Error in adding cycle(s) to potential: #coefficients != #cycles\nFailed to add potential")
            else:
                for iec, ec in enumerate(edge_cycle):
                    c = coef[iec]
                    cycle=[]
                    for i,vi in enumerate(ec):
                        vj = ec[(i+1)%len(ec)]
                        try:
                            cycle.append(self.Q1.index([vi,vj]))
                        except:
                            print("\nError in add_term_to_potential: there is no edge with endpoints (%d, %d). Ignoring term\n"%(vi,vj))
                    cycle = cycleOrder(tuple(cycle))
                    self.potential[cycle] = coef[iec]


    def __repr__(self):
        return "quiver with potential:\nedges: " \
               +repr(self.Q1) \
               +"\npotential:"+self.print_potential()


    def add_edges(self, edge):
        ctr = len(self.Q1)
        for e in edge:
            self.Q1.append(e)
            if e[0] == e[1]:
                self.loops_at[e[0]].append((ctr, [e[0],e[0]]))
                self.can_mutate[e[0]] = False
            else:
                self.arrows_with_head[e[1]].append((ctr,e))
                self.arrows_with_tail[e[0]].append((ctr,e))
            ctr += 1
        self.incidence_matrix = matrixFromEdges(self.Q1)


    def remove_edges(self, edge):
        # create list of edges to keep
        keep = [i for i in range(len(self.Q1))]
        for e in edge:
            keep[e] = -(e+1)

        self.Q1 = [self.Q1[i] for i in keep if i >= 0]

        k = [x for x in keep if x >= 0]
        new_edge_indices = dict(zip(k, range(len(k))))
        for v in self.Q0:
            self.arrows_with_head[v] = [(new_edge_indices[i], e) for (i,e) in self.arrows_with_head[v] if keep[i]>=0]
            self.arrows_with_tail[v] = [(new_edge_indices[i], e) for (i,e) in self.arrows_with_tail[v] if keep[i]>=0]
            self.loops_at[v] = [(new_edge_indices[i], self.Q1[new_edge_indices[i]]) for (i,e) in self.loops_at[v] if keep[i]>0]
        self.incidence_matrix = matrixFromEdges(self.Q1)

        self.can_mutate = [True if len(self.loops_at[x]) < 1 else False for x in self.Q0]
        #update the potential 
        p = {}
        for term, coef in self.potential.items():
            try:
                new_term = tuple([new_edge_indices[x] for x in term])
                p[new_term] = coef
            except:
                print("\nError in remove_edges: one of the edges to remove shows up in potential. Dropping that term\n")
        self.potential = p


    def mutate(self, v, warnings=True):
        # make a copy of the original quiver to mutate
        QP = copy.deepcopy(self)
    
        # first check if there's a loop at v. 
        if not self.can_mutate[v]:#len(self.loops_at[v]) > 0:
            if warnings:
                print("Error in mutate routine: there's a loop at vertex %d. returning"%v)
            return QP

        # reverse all edges incident to v:
        for (i, e) in self.arrows_with_tail[v] + self.arrows_with_head[v]:
            QP.Q1[i] = [e[1],e[0]]

        # update connectivity info for quiver
        QP.arrows_with_head[v] = list([(i, QP.Q1[i]) for (i, e) in self.arrows_with_tail[v]])
        QP.arrows_with_tail[v] = list([(i, QP.Q1[i]) for (i, e) in self.arrows_with_head[v]])

        delta = {}
        shortcuts = {}
        new_edges = []
        edgectr = len(QP.Q1)
        # add 'shortcuts' i->j for all 2-paths(which are then reversed): i->v->j 
        # and keep track of the 3-cycles that are produced (saved in list delta)
        for e1i, e1 in self.arrows_with_head[v]:
            i = e1[0]
            for e2i, e2 in self.arrows_with_tail[v]:
                j = e2[1]

                # add the shortcut edge
                new_edges.append([i,j])

                # keep track of 3-cycle introduced by adding this edge
                delta[(e2i, e1i, edgectr)] = 1
                shortcuts[tuple([e1i,e2i])] = edgectr
                edgectr += 1
        QP.add_edges(new_edges)

        # update the potential
        wprime = {}
        # first replace each 2-path i->v->j with its shortcut i->j
        for monoid, coef in self.potential.items():
            m = []
            for i, x1 in enumerate(monoid):
                if (x1, self.Q1[x1]) in self.arrows_with_head[v]:

                    x2 = monoid[(i+1)%len(monoid)]
                    pair = tuple([x1, x2])

                    if pair in shortcuts.keys():
                        m.append(shortcuts[pair])
                    else:
                        m.append(x1)
                        m.append(x2)
                else:
                    x2 = monoid[(i+len(monoid)-1)%len(monoid)]
                    if (x2, self.Q1[x2]) not in self.arrows_with_head[v]:
                        m.append(x1)
            wprime[tuple(m)] = coef

        # add the set of 3-cycles introduced with the shortcuts
        QP.potential = {**wprime, **delta}

        return reduce_QP(QP, warnings=warnings)


    def draw(self, time=1, **kwargs):
        g = nx.DiGraph()
        for x in self.Q0:
            g.add_node(str(x))
        for e in self.Q1:
            if e[0] != e[1]:
                g.add_edge(str(e[0]), str(e[1]))
        pos = nx.spring_layout(g)
        if self.positions is not None:
            pos = {str(i):v for i,v in enumerate(self.positions)}

        nx.draw_networkx(g, pos, with_labels=True, connectionstyle='arc3, rad=0.1', node_size=200, **kwargs)
        for v in self.Q0:
            for i in range(len(self.loops_at[v])):
                new_edge = [(str(v), str(v))]
                g.add_edges_from(new_edge)
                nx.draw_networkx_edges(g, pos, edgelist=new_edge, connectionstyle='arc3, rad=0.1', node_size=(i+2)*200)

        txt = "W = " + self.print_potential()
        plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center')
        plt.draw()
        plt.pause(time)
        plt.clf()


    def mutate_in_sequence(self, vertex_sequence, draw=True):
        other = copy.deepcopy(self)
        if len(vertex_sequence) > 0:
            for v in vertex_sequence:
                QP = other.mutate(v)
                vertex_colors = ['b' if i != v else 'r' for i in other.Q0]
                if draw:
                    kw = {'node_color': vertex_colors}
                    other.draw(time=.75, **kw)
                    QP.draw(time=.75, **kw)
                other=QP
        return other




class Resolution():

    def __init__(self, e, vertex_positions=None):
        self.edges = e
        triangulation = make_triangulation(self.edges)
        self.triangles = triangulation[0]
        self.edge_to_triangle = triangulation[1]
        self.vertex_positions = vertex_positions


    def update(self, R):
        self.edges = copy.deepcopy(R.edges)
        self.triangles = copy.deepcopy(R.triangles)
        self.edge_to_triangle = copy.deepcopy(R.edge_to_triangle)
        self.vertex_positions = copy.deepcopy(R.vertex_positions)


    def QP(self):
        edges = [(t[i%3], t[(i+j)%3]) for j in [-1,1] for i in range(1, 4) for t in self.triangles]
        for i, e in enumerate(self.edges):
            if self.curve_type(e) == (-2,0):
                edges.append((i,i))
            if self.curve_type(e) == (-3,0):
                edges.append((i,i))
                edges.append((i,i))

        p = potential(edges, self.triangles)

        positions = None
        if self.vertex_positions is not None:
            positions = [[0.5*(self.vertex_positions[e[0]][0]
                              +self.vertex_positions[e[1]][0]), 
                          0.5*(self.vertex_positions[e[0]][1]
                              +self.vertex_positions[e[1]][1])] for e in self.edges]
        return QuiverWithPotential(edges, potential=p, positions=positions)


    def draw(self, time=0.5, **kwargs):
        # create networkx graph for drawing purposes
        g = nx.Graph()
        g.add_edges_from(self.edges)

        # option for coloring specific edges
        for key in kwargs.keys():
            if 'edge' in key:
                # have to manually match entries of edge_colors since nx reorders edges when it creates a graph
                for e in g.edges:
                    if list(e) in self.edges:
                        i = self.edges.index(list(e))
                    else:
                        i = self.edges.index([e[1],e[0]])
                    g.edges[e][key] = kwargs[key][i]
                kwargs[key] = [x for x in nx.get_edge_attributes(g,key).values()]

        # option for drawing with specified vertex positions
        if self.vertex_positions is not None:
            nx.draw_networkx(g, self.vertex_positions, **kwargs)
        else:
            nx.draw_networkx(g, **kwargs)

        plt.draw()
        plt.pause(time)
        plt.clf()


    def can_flop(self, e):
        """
        check if an edge e (passed as either the tuple representation or its index in the resolution) can be flopped
        """
        return self.curve_type(e) == (-1,-1)


    def curve_type(self, e):
        if isinstance(e, list):
            ei = self.edges.index(e)
        else:
            ei, e = e, self.edges[e]

        e2t = sorted(self.edge_to_triangle[ei])
        if len(e2t) < 2: # check if its an edge on the boundary
            return (0,0)

        # unpack the vertices on the edge and those it would be flopped to 
        vertices = set([tuple(self.vertex_positions[i]) for t in e2t for e in self.triangles[t] for i in self.edges[e]])
        # vertices of existing edge
        e_verts = [tuple(self.vertex_positions[e[0]]), tuple(self.vertex_positions[e[1]])]

        # vertices and direction vector of flopped edge
        o_verts = [x for x in vertices if x not in e_verts]
        alt_e = [o_verts[1][0] - o_verts[0][0], o_verts[1][1] - o_verts[0][1]]

        """
         e_verts[0]     o_verts[1]
         *-------------------*
         | \              /  |
         |   \     alt_e.    |
         |     \      /      |
         |       \  .        |
         |         \         |
         |t1v1   .   \       |
         |     /      e\     |
         |   .           \   |
         | /    t2v1       \ |
         *___________________*
         o_verts[0]      e_verts[1]


        Note that if e can be flopped(i.e. if the two triangles form a
        convex 4-sided object), then the angles:
          * angle1 (between t2v1 and alt_e)
          * angle2 (between alt_e and t1v1)
        will be both positive (if e_verts is oriented as in the
        picture above), or both negative (if e_verts is reversed).
        """

        t1v1 = [e_verts[0][0] - o_verts[0][0], e_verts[0][1] - o_verts[0][1]]
        t2v1 = [e_verts[1][0] - o_verts[0][0], e_verts[1][1] - o_verts[0][1]]

        # affine transformation so that arctan will avoid branch cuts
        # when computing the angles between alt_e and the two sides
        # t1v1 and t1v2.
        v1 = rotated_vector(basis=t1v1, hyp=alt_e)
        v2 = rotated_vector(hyp=t2v1, basis=alt_e)
        angle1 = np.arctan2(v1[1], v1[0])
        angle2 = np.arctan2(v2[1], v2[0])

        tol = 1.0e-6
        if np.sign(angle1) == np.sign(angle2):
            # check that neither angle is zero
            if (abs(angle1) > tol) and (abs(angle2) > tol):
                return (-1,-1)
            return (-2,0)
        return (-3,0)


    def flop_in_sequence(self, edge_sequence, draw=True):
        for edge in edge_sequence:
            R = self.flop(edge)
            edge_colors = ['b' if i != edge else 'r' for i in range(len(self.edges))]
            if draw:
                kw = {'edge_color': edge_colors}
                self.draw(time=.5, **kw)
                R.draw(time=.75, **kw)
            self.update(R)
        return R


    def __repr__(self):
        return "resolution:\n\t%d triangles\n\t"%(len(self.triangles))+repr(self.triangles)


    def flop(self, e):
        # make a copy of self to return: 
        toRet = copy.deepcopy(self)

        # check that this edge can be flopped 
        # (i.e. it's not on a boundary, and the quadrilateral
        # formed by the edge's neighboring triangles is convex)
        if not self.can_flop(e):
            return toRet

        # find the two triangles that share edge e
        nbrs = self.edge_to_triangle[e]
        nbr1,nbr2 = nbrs[:2]

        # get the points p1, p2 'opposite' to edge e
        #    p1-------e[1]
        #    |       /  |
        #    | nbr1 /   |
        #    |     /e   |
        #    |    /     |
        #    |   / nbr2 |
        #    |  /       |
        #   e[0]-------p2
        p1 = [y for x in self.triangles[nbr1] for y in self.edges[x] if y not in self.edges[e]][0]
        p2 = [y for x in self.triangles[nbr2] for y in self.edges[x] if y not in self.edges[e]][0]

        # find the edge indices for e[0]--p2 and e[1]--p1 and other 2 as well.
        e1p1 = [y for y in self.triangles[nbr1] if ((self.edges[e][1] in self.edges[y]) and (y != e))][0]
        e0p1 = [y for y in self.triangles[nbr1] if ((self.edges[e][0] in self.edges[y]) and (y != e))][0]
        e0p2 = [y for y in self.triangles[nbr2] if ((self.edges[e][0] in self.edges[y]) and (y != e))][0]
        e1p2 = [y for y in self.triangles[nbr2] if ((self.edges[e][1] in self.edges[y]) and (y != e))][0]

        # update the triangle edges
        toRet.triangles[nbr1] = [e0p2,e,e0p1]
        toRet.triangles[nbr2] = [e1p2,e,e1p1]
        toRet.edges[e] = [p2, p1]

        # replace the connectivity for edges to triangles
        toRet.edge_to_triangle[e1p1] = [y if y != nbr1 else nbr2 for y in self.edge_to_triangle[e1p1]]
        toRet.edge_to_triangle[e0p2] = [y if y != nbr2 else nbr1 for y in self.edge_to_triangle[e0p2]]

        return toRet


def potential(edges, triangles):
    # build triples. Note: since the triangles' edges aren't necessarily ordered, we need to impose order here.
    triples = [[6*ti + x for x in i] for i in [[1,3,5],[0,2,4]] for ti in range(len(triangles))]
    p = {}
    for t in triples:
        t1,h1 = edges[t[0]]
        t2,h2 = edges[t[1]]
        if h1 == t2:
            orderedTriple = tuple(t)
        else:
            orderedTriple = (t[1],t[0],t[2])
        p[cycleOrder(orderedTriple)] = 1
    return p


def make_triangulation(edges):
    """
        this routine takes a list of edges where each edge is in the format [v1, v2] 
        and creates an edge-to-triangle and triangle-to-edge connectivity structure
        for all triangles that can be formed based on vertex identification. 
    """

    numVerts = len(set(flattenList(edges)))
    nbrs = [set() for x in range(numVerts)] # create a set of edges incident to each vertex (empty at first)
    for i, e in enumerate(edges): # populate the list of edges incident to each vertex
        nbrs[e[0]].add(i)
        nbrs[e[1]].add(i)

    triangles = []
    edge_to_triangle = [[] for i in edges]

    triples = set()
    triangle_ctr = 0
    for edge1, e in enumerate(edges):
        # extract the vertices at ends of edge1
        h,t=e
        # find all of the edges connected to these endpoints
        h_nbrs = nbrs[h]
        t_nbrs = nbrs[t]

        # loop over all edges that connect to one endpoint of edge1
        for edge2 in t_nbrs:
            if edge2 != edge1:
                p0, p1 = edges[edge2] # get the endpoints of neighboring edge
                p = p1 if (p0 == t) else p0

                # if the neighboring edge connects to one of the neighbors
                # of the second endpoint, then we have a triangle.
                e3 = nbrs[p].intersection(h_nbrs)
                if len(e3) > 0:
                    edge3 = e3.pop()

                    # check that the edge triple hasn't already been added
                    st = ",".join([str(x) for x in sorted([edge1,edge2,edge3])])
                    if st not in triples:
                        triangles.append([edge1, edge2, edge3])

                        edge_to_triangle[edge1] += [triangle_ctr]
                        edge_to_triangle[edge2] += [triangle_ctr]
                        edge_to_triangle[edge3] += [triangle_ctr]
                        triples.add(st)
                        triangle_ctr += 1

    return triangles, edge_to_triangle


def path_derivative(potential, arrow):
    """
    take the path derivative of each term of the potential w.r.t. the arrow provided
    """
    to_return = {}
    for term, coef in potential.items():
        if arrow in term:
            i = term.index(arrow)
            new_term = term[i+1:]+term[:i]
            to_return[new_term] = coef
    return to_return


def edgesFromMatrix(mat):
    """
    creates a list of edges from an incidence matrix
    """
    return [[r.index(-1), r.index(1)] for r in np.matrix(mat).transpose().tolist()]


def matrixFromEdges(edges, oriented=True):
    """
    create an incidence matrix from a list of edges
    """
    nv = len(set(np.array(edges).flatten()))
    tail = -1 if oriented else 1
    m = np.matrix([[tail if x == e[0] else 1 if x == e[1] else 0 for x in range(nv)] for e in edges]).transpose()
    for i,e in enumerate(edges):
        if e[0] == e[1]:
            m[e[0],e] = 2
    return m


def reduce_QP(QP, warnings=True):
    """
        this routine takes a QP and "reduces" it by removing all 2-cycles that 
        show up in the potential W (as long as each term in the 2-cycle can be 
        removed, which depends on W), and also removing the associated edges in QP.Q1
    """
    # compute the partial derivative of QP's potential for every edge
    partials = [path_derivative(QP.potential, a) for a in range(len(QP.Q1))]

    # create a lookup of the replacements associated to each partial derivative
    # this is based on the fact that equivalent terms in the jacobian algebra can 
    # be found by factoring the terms of dW/da for each arrow a. 
    reduce_dict = {}
    for pd in partials:
        for k,v in pd.items():
            reduce_dict[k] = [(k1,-v1/v) for k1,v1 in pd.items() if k1 != k]

    # find out which edges show up in quadratic terms in the potential. 
    edges_to_remove = sorted([x for term in QP.potential.keys() for x in term if (len(list(term)) == 2)])
    for e in range(len(QP.Q1)):
        if e not in edges_to_remove:
            reduce_dict[tuple([e])] = [(e,1)]

    # find out which edges are equivalent to 0
    zero_terms = [k for (k,v) in reduce_dict.items() if len(v) < 1]

    # now create a lookup table of replacements for each edge to be removed
    # and each replacement term only contains edges that are NOT in edges_to_remove+zero_terms
    for e in edges_to_remove:

        # unpack the list of monoids that sum up to something equivalent to e in the jacobian algebra
        terms_for_e = reduce_dict[tuple([e])]
        if len(terms_for_e) > 0:

            found_replacement=False; ctr = 0 # stopping criterion
            while not found_replacement and ctr < len(edges_to_remove):
                found_replacement=True; ctr += 1 # update stopping criterion

                # placeholder for terms to be added
                alt_terms_for_e = []
                for (term, coef) in terms_for_e:
                    alt_term = [(term, coef)]

                    # check if any of the terms in e's replacement terms 
                    # also contains one of the edges_to_remove
                    if any([x in edges_to_remove for x in term]):
                        found_replacement=False

                        # if so, then we need to replace that term. 
                        new_term = [((),1)]
                        for tti, tt in enumerate(term):
                            if tt in edges_to_remove:
                                new_term = [(a+b, c*d) for (a,c) in new_term \
                                        for (b,d) in reduce_dict[tuple([tt])] \
                                        if len(reduce_dict[tuple([tt])]) > 0]
                            else:
                                new_term = [(a+tuple([tt]), c) for (a,c) in new_term]
                        alt_term = [(a, c*coef) for (a,c) in new_term]
                    alt_terms_for_e.extend(alt_term)
                terms_for_e = alt_terms_for_e

            reduce_dict[tuple([e])] = terms_for_e
            if not found_replacement and warnings:
                print("problem in reducing QP: could not properly remove edge %d"%e)

    # now create a global lookup: each edge has a "replacement" 
    # (identity or set of terms, depending on whether or not it's in edges_to_remove)
    new_reduce_dict = {e: [(tuple([e]),1)] for e in range(len(QP.Q1))}
    for e in edges_to_remove:
        new_reduce_dict[e] = reduce_dict[tuple([e])]

    # keep track of which edges are equivalent to 0. 
    zero_terms = []
    for k,v in new_reduce_dict.items():
        if len(v) < 1:
            zero_terms.append(k)

    # now update the potential by replacing all of the terms in edges_to_remove
    Wprime = {}
    for term, coef in QP.potential.items():
        if not any([tuple([x]) in zero_terms for x in term]):
            new_terms = [((), 1)]
            for tt in term:
                new_terms = [(a+b, c*d) for (a,c) in new_terms \
                        for (b,d) in new_reduce_dict[tt] \
                        if len(new_reduce_dict[tt]) > 0]
            new_terms = [(a, b*coef) for (a,b) in new_terms]
            for (y,z) in new_terms:
                if cycleOrder(y) in Wprime:
                    Wprime[cycleOrder(y)] += z
                else:
                    Wprime[cycleOrder(y)] = z

    QP.potential = {k:v for k, v in Wprime.items() if v != 0}
    QP.remove_edges(sorted(edges_to_remove))
    return QP


def cycleOrder(cycle):
    cycle = tuple(cycle)
    minidx = cycle.index(min(cycle))
    return tuple(cycle[minidx:]+cycle[:minidx])


def rotated_vector(basis=[1,0], hyp=[1,1]):
    """
    Inputs:
    This routine takes two vectors basis and hyp, both
    of which are assumed to have tails at the origin,
    so only the head is passed in each case.

    Outputs:
    The computes the affine transformation necessary to map
    The input "basis" to <1,0>, and returns the image of "hyp"
    under that transformation.
    """
    bb = np.array([basis[0], basis[1]])
    hh = np.matrix([hyp[0], hyp[1]]).transpose()

    # create matrix that translates basis to the line segment [0,1]
    if bb[0] != 0:
        matrix_of_transformation = np.matrix([[-bb[1]+1/bb[0],bb[0]],[-bb[1],bb[0]]])
    else:
        matrix_of_transformation = np.matrix([[-bb[1],bb[0]+1/bb[1]],[-bb[1],bb[0]]])

    return (matrix_of_transformation*hh).transpose().tolist()[0]


def current_QP(QP):
    """ return a string representation of the QP (used to check
    uniqueness when calculating all possible mutations of a 
    given QP. More efficient than built-in equality ftn) """
    sorted_items = sorted(["%d_%d_"%(e[0],e[1]) for e in QP.Q1]) + ["|"] \
                 + sorted([".%d_"%int(v)+"_".join(["%d"%tk for tk in k]) for k,v in QP.potential.items()])
    return "".join(sorted_items)


def calculate_all_mutations_from_vertex(Q, v=0, already_met=set(), saved = []):
    """calculates all of the unique Quivers obtained from the mutation of 
    input quiver Q at the vertex v"""
    met_the_end = True
    mutation = Q.mutate(v, warnings=False)
    cQ = current_QP(Q)
    cm = current_QP(mutation)

    if cQ[:cQ.index("|")] not in already_met:
        saved.append(cQ)
    if cm[:cm.index("|")] not in already_met:
        saved.append(cm)

    already_met.add(cQ[:cQ.index("|")])

    # if we've already met this mutation or we can't mutate
    if (cm[:cm.index("|")] in already_met) or (len(Q.loops_at[v]) > 0):
        return True, already_met, saved

    # otherwise try mutating new QP at every other vertex
    already_met.add(cm[:cm.index("|")])
    for other_v in mutation.Q0:
        met_an_end, already_met, saved = calculate_all_mutations_from_vertex(mutation, other_v, already_met)
        if not met_an_end:
            met_the_end = False

    return met_the_end, already_met, saved



def get_all_mutations_from_quiver(QP):
    """This routine calculates all possible quivers that can be obtained 
    by mutating any set of vertices of QP.
        inputs: QP = Quiver with Potential
        outputs: list of QPs (all possible unique quivers)
    """
    all_ms = set()
    for v in QP.Q0:
        tf, sms, ms = calculate_all_mutations_from_vertex(QP, v)
        all_ms = all_ms.union(set(ms))
    ms = list(all_ms)

    for m in ms:
        e, p = m.split("|")
        edges = e.split("_")
        p = p.split(".")

        edges = [[int(edges[2*i+j]) for j in [0,1]] for i in range(int(len(edges)/2))]

        pot = {}
        for v in p:
            if len(v) > 0:
                vs = v.split("_")
                pot[(int(x) for x in vs[1:])] = int(v.split("_")[0])

        yield QuiverWithPotential(edges, pot, QP.positions)


def calc_mutations_by_index_list(Q, v=0, already_met=set(), current_path = [], saved_paths = []):
    """calculates all of the unique quivers obtained from the mutation of 
    input quiver Q at the vertex v but returns the (non unique) sequence
    of vertices to mutate in order to obtain each unique quiver"""
    met_the_end = True
    mutation = Q.mutate(v, warnings=False)
    cQ = current_QP(Q)
    cm = current_QP(mutation)

    already_met.add(cQ[:cQ.index("|")])
    saved_paths.append(list(current_path))

    # if we've already met this mutation or we can't mutate
    if (cm[:cm.index("|")] in already_met) or (len(Q.loops_at[v]) > 0):
        return True, already_met, current_path, saved_paths

    current_path = current_path + [v]
    # otherwise try mutating new QP at every other vertex
    already_met.add(cm[:cm.index("|")])
    for other_v in mutation.Q0:
        this_path = list(current_path)
        met_an_end, already_met, this_path, saved_paths = calc_mutations_by_index_list(mutation, other_v, already_met, this_path, saved_paths)

        if not met_an_end:
            met_the_end = False

    return met_the_end, already_met, current_path, saved_paths


def all_mutation_sequences_for_quiver(Q):
    """This routine calculates all possible quivers that can be obtained 
    by mutating any set of vertices of QP.
        inputs: QP = Quiver with Potential
        outputs: list of (ordered) sequences of vertices to mutate in order to obtain all possible unique quivers
    """
    all_ms = []
    for v in Q.Q0:
        tf, sms, cms, ms = calc_mutations_by_index_list(Q, v)
        all_ms = all_ms + ms

    reduced_ms = []
    already_met = set()
    for m in all_ms:
        qp = Q.mutate_in_sequence(m, draw=False)
        cm = current_QP(qp)
        if cm[:cm.index("|")] not in already_met:
            already_met.add(cm[:cm.index("|")])
            reduced_ms.append(m)
    return reduced_ms
