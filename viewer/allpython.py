import numpy as np
import copy
import networkx as nx
import matplotlib.pyplot as plt
import json

def flattenList(l):
    if type(l) is list or type(l) is tuple:
        return [x for y in l for x in flattenList(y)]
    return [l]


class QuiverWithPotential():

    def __init__(self, graph_obj, potential=None, positions=None, frozen_nodes=[]):
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
        self.frozen_nodes = frozen_nodes

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
        for n in frozen_nodes:
            self.can_mutate[n] = False


    def print_potential(self):
        return " + ".join(["%.3lf X_"%v+".".join([str(self.Q1[x][0]) for x in k]) for k,v in self.potential.items()])


    def add_term_to_potential(self, edge_cycle, coef=1, input_format="vertex_order"):
        if not isinstance(edge_cycle[0], tuple):
            edge_cycle = [edge_cycle]

        if coef == 1:
            coef = [1 for e in edge_cycle]
        if len(coef) != len(edge_cycle):
            print("Error in adding cycle(s) to potential: #coefficients != #cycles\nFailed to add potential")

        cycles=[]
        if input_format == "vertex_order":
            for ec in edge_cycle:
                cycle=[]
                for i,vi in enumerate(ec):
                    vj = ec[(i+1)%len(ec)]
                    try:
                        cycle.append(self.Q1.index([vi,vj]))
                    except:
                        print("\nError in add_term_to_potential: there is no edge with endpoints (%d, %d). Ignoring term\n"%(vi,vj))
                cycle = cycleOrder(tuple(cycle))
                cycles.append(cycle)
            edge_cycle=cycles

        for i, e in enumerate(edge_cycle):
            self.potential[cycleOrder(e)] = coef[i]


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
                self.arrows_with_head[e[1]].append((ctr,[e[0],e[1]]))
                self.arrows_with_tail[e[0]].append((ctr,[e[0],e[1]]))
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
        for n in self.frozen_nodes:
            self.can_mutate[n] = False

        #update the potential 
        p = {}
        for term, coef in self.potential.items():
            try:
                new_term = tuple([new_edge_indices[x] for x in term])
                p[new_term] = coef
            except:
                print("\nError in remove_edges: one of the edges to remove shows up in potential. Dropping that term\n")
        self.potential = p


    def remove_nonloop_multiedges(self):
        # first group edges by which vertices they adjoin
        replacing = [[i for i,v in self.arrows_with_head[e[1]] \
                        if v[0] == e[0]] \
                      for e in self.Q1 if e[0] != e[1]]
        # create a lookup: replace_with[y] = x means that y gets replaced with x
        replace_with = {y: min(x) for x in replacing for y in x}
        for ei in range(len(self.Q1)):
            if ei not in replace_with.keys():
                replace_with[ei] = ei

        # now update the potential
        new_potential = {}
        for k,v in self.potential.items():
            new_term = tuple([replace_with[x] for x in k])
            new_potential[new_term] = v
        self.potential = new_potential
        remove_edges = sorted([y for y,r in replace_with.items() if ((y!=r) and (self.Q1[y][0] != self.Q1[y][1]))])

        # finally, remove the edges that have been replaced
        self.remove_edges(remove_edges)


    def mutate(self, v, warnings=True):
        # make a copy of the original quiver to mutate
        QP = copy.deepcopy(self)
    
        # first check if there's a loop at v. 
        if not self.can_mutate[v]:
            if warnings:
                print("Error in mutate routine: cannot mutate at vertex %d. returning"%v)
            return QP

        saved_edges = [[e[0],e[1]] for e in QP.Q1]
        # reverse all edges incident to v:
        for (i, e) in self.arrows_with_tail[v] + self.arrows_with_head[v]:
            saved_edges[i] = [e[1],e[0]]

        delta = {}
        shortcuts = {}
        new_edges = []
        edgectr = len(self.Q1)
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
        saved_edges.extend(new_edges)

        # update the potential
        wprime = {}
        # first replace each 2-path i->v->j with its shortcut i->j
        for monoid, coef in self.potential.items():
            m = []
            for i, x1 in enumerate(monoid):
                if saved_edges[x1][0] == v:

                    x2 = monoid[(i+1)%len(monoid)]
                    pair = tuple([x1, x2])

                    if pair in shortcuts.keys():
                        m.append(shortcuts[pair])
                    else:
                        m.append(x1)
                        m.append(x2)
                else:
                    x2 = monoid[(i+len(monoid)-1)%len(monoid)]
                    if saved_edges[x2][0] != v:
                        m.append(x1)
            wprime[tuple(m)] = coef

        # add the set of 3-cycles introduced with the shortcuts
        saved_potential = {**wprime, **delta}
        saved_positions = copy.deepcopy(self.positions)
        QP = QuiverWithPotential(saved_edges, potential=saved_potential, positions=saved_positions, frozen_nodes = [x for x in self.frozen_nodes])
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


    def mutate_in_sequence(self, vertex_sequence, draw=True, warnings=True):
        other = copy.deepcopy(self)
        if len(vertex_sequence) > 0:
            for v in vertex_sequence:
                QP = other.mutate(v, warnings=warnings)
                vertex_colors = ['b' if i != v else 'r' for i in other.Q0]
                if draw:
                    kw = {'node_color': vertex_colors}
                    other.draw(time=.75, **kw)
                    QP.draw(time=.75, **kw)
                other=QP
        return other

    def toJSON(self, filename=None):
        ns = [{"id":str(ni), \
               "label":str(ni), \
               "x":100*xy[0], \
               "y":100*xy[1]} \
               for ni, xy in enumerate(self.positions)
               ]
        es = [{"id":str(ei), \
               "title": "edge %d"%ei, \
               "from": str(e[0]), \
               "to": str(e[1]), \
               "arrows": "to"} \
               for ei, e in enumerate(self.Q1)
               ]
        fn = [{"id":str(x)} for x in self.frozen_nodes]
        pt = [{"id":",".join(k), "coef":str(v)} for k,v in self.potential.items()]
        obj = json.dumps({"nodes": ns, "edges":es, "frozenNodes":fn, "potential":pt}, indent=4)
        if filename is not None:
            with open(filename, "w") as f:
                f.write(obj)
        return obj



def make_triangulation(edges, extremal_edges=None):
    """
        this routine takes a list of edges where each edge is in the format [v1, v2] 
        and creates an edge-to-triangle and triangle-to-edge connectivity structure
        for all triangles that can be formed based on vertex identification. 
    """

    numVerts = max(flattenList(edges)) + 1
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

    # remove the triangles that are extra, and get included because they are subdivisions. 
    # That is, triangles like the boundary of the following:
    """
                        *
                      / | \
                     /  |  \
                    /  .*.  \
                   / ./   \. \
                  /./       \.\
                 *-------------*
    """
    if len(triangles) > 1:
        if extremal_edges is not None:
            if len(extremal_edges) == 3:
                extremal_triangle = \
                        set(edge_to_triangle[extremal_edges[0]]).intersection(
                        set(edge_to_triangle[extremal_edges[1]])).intersection(
                        set(edge_to_triangle[extremal_edges[2]])).pop()
                triangles = [t for ti, t in enumerate(triangles) if ti != extremal_triangle]
                edge_to_triangle = [[t - int(t > extremal_triangle) for t in e2t if t != extremal_triangle] for e2t in edge_to_triangle]

        ts_to_keep = []
        for ti, t in enumerate(triangles):
            if all([((len(edge_to_triangle[e]) < 2) or (len(edge_to_triangle[e]) > 2)) for e in t]):
                pass
            else:
                ts_to_keep.append(ti)
        triangles = [triangles[i] for i in ts_to_keep]
        edge_to_triangle = [[ts_to_keep.index(t) for t in e2t if t in ts_to_keep] for e2t in edge_to_triangle]

    return triangles, edge_to_triangle


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


def find_dependencies(current_item, met, lookups):
    if current_item in met:
        return True, met

    met = met + [current_item]
    for item in lookups[current_item]:
        return find_dependencies(item, met, lookups)
    return False, met



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
    edges_to_remove = sorted(list(set([x for term in QP.potential.keys() for x in term if (len(list(term)) == 2)])))

    # check if there is a circular dependency between any of the quadratic terms 
    # (this makes it so that one edge at least cannot be deleted in the reduction process)
    if len(edges_to_remove)%2 > 0:
        lookups = {x:[z \
                for y in reduce_dict[tuple([x])] for z in y[0] \
                if (z in edges_to_remove) \
                ] for x in edges_to_remove
                }
        loops = set()
        for e in edges_to_remove:
            in_a_loop, loop = find_dependencies(e, [], lookups)
            if in_a_loop:
                loops.add(tuple(sorted(loop)))
        if len(loops) > 0:
            quad_edges_to_keep = [x[0] for x in loops]
            edges_to_remove = [x for x in edges_to_remove if x not in quad_edges_to_keep]

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
        this_path = copy.deepcopy(current_path)
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
        qp = Q.mutate_in_sequence(m, draw=False, warnings=False)
        cm = current_QP(qp)
        if cm[:cm.index("|")] not in already_met:
            already_met.add(cm[:cm.index("|")])
            reduced_ms.append(m)
    return reduced_ms


def lattice_points_nonunit(r=6,a=1,b=2,c=3):
    # Lattice points in the junior simplex lie in the plane x+y+z=1 and are of the form 
    # (i/r,j/r,k/r) where i+j+k=r. 
    # This routine returns the triples (i,j,k) instead of scaling by 1/r. 
    points = [(r,0,0), (0,r,0), (0,0,r)]
    for i in range(r):
        ai = (a*i)%r
        bi = (b*i)%r
        ci = (c*i)%r
        if (ai + bi + ci) == r:
            points.append((ai,bi,ci))

    return list(set(points))

def convex_hull(planar_points):
    # performs a QuickHull search to find which points in the set of 
    # planar_points form the convex hull. 
    # This routine assumes that the set of planar_points are 
    # points in any dimension, but which lie in some plane, and that the 
    # first two points in the set planar_points form an edge of the 
    # convex hull. The general version defines the initial starting points 
    # p0 and p1 as any two points whose distance is maximal among pairs of points

    def findHull(pts, a, b):
        # line segment from a to b and its length
        AB = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
        rays = [np.array(p) - np.array(a) for p in pts]

        # find all the points on one side of segment ab
        pts = [np.array(p) for i, p in enumerate(pts) if np.all(np.cross(AB, rays[i]) >= 0)]
        rays = [np.array(p) - np.array(a) for p in pts]

        if len(pts) > 0:
            # calculate length of segment ab
            ABL = np.dot(AB, AB)

            # get length of vector from each pt to line segment AB
            distances = [r - (np.dot(AB,r)/ABL)*AB for r in rays]
            dLs = [np.dot(d,d) for d in distances]
            maxPt = dLs.index(max(dLs))

            dpts = pts[:maxPt]+pts[maxPt+1:]
            return findHull(dpts, a, pts[maxPt]) + findHull(dpts, pts[maxPt], b)
        else:
            return [tuple(a), tuple(b)]

    p0,p1 = planar_points[0],planar_points[1]
    points = list(set(findHull(planar_points[2:], p0, p1) + findHull(planar_points[2:], p1, p0)))
    planar_points = [tuple(p) for p in planar_points]
    return sorted([planar_points.index(p) for p in points])


def ordered_rays(L):
    # this routine assumes that there exists a Jung–Hirzebruch 
    # continued fraction relation among the rays in the list L
    # and returns the order and 'strength' of each ray in the 
    # form of a list with entries of the form (indexi, ai)
    # to denote the index of L[i] in the ordered list, and its strength ai
    # NOTE: it also assumes that the first and last entry in L 
    # corresponds each to a ray at one 'end'

    L = [np.array(l) for l in L]
    l0,lk = L[0],L[-1]
    iLs = L[1:-1]
    if len(iLs) < 1:
        return [(0, 0), (1, 0)]
    else:
        angles = np.argsort([veclen(np.cross(l, l0)) for l in iLs])
        oLs = [l0] + [iLs[x] for x in angles] + [lk]
        toRet = [(0,0)] + [(angles[i-1]+1, (oLs[i-1][np.flatnonzero(oLs[i])[0]] + oLs[i+1][np.flatnonzero(oLs[i])[0]]) / oLs[i][np.flatnonzero(oLs[i])[0]]) for i in range(1, len(oLs)-1)] + [(len(L)-1,0)]
        return toRet


def is_collinear(ray, pt, tol=1.0e-6):
    # this routine checks if a point lies along the line 
    # defined by the two n-dimensional points ray = [point1, point2]

    # first check if the point is 'equal' to one of the points defining the ray
    if np.allclose(np.array(ray[0]), np.array(pt)) or np.allclose(np.array(ray[1]), np.array(pt)):
        return True

    # otherwise, we'll look at the component-wise slopes 
    numerator = np.array(pt) - np.array(ray[0])
    denominator = np.array(ray[1]) - np.array(ray[0])

    if set(numerator.nonzero()[0]) != set(denominator.nonzero()[0]):
        return False

    nz = np.nonzero(denominator)
    frac = numerator[nz] / denominator[nz]
    return frac.max() - frac.min() < tol


def intersects(seg1, seg2, pts):
    # check if two segments intersect at a point interior to both 
    # (i.e. excluding intersection at endpoints and intersection
    # of the extended lines passing through the segments)
    a1 = np.array(pts[seg1[0]])
    a2 = np.array(pts[seg2[0]])
    b1 = np.array(pts[seg1[1]]) - a1
    b2 = np.array(pts[seg2[1]]) - a2

    if np.count_nonzero(np.cross(b1,b2)) >= len(b1):
        t1 = np.cross(a2-a1, b2) / np.cross(b1,b2)
        t2 = np.cross(a1-a2, b1) / np.cross(b2,b1)

        if np.count_nonzero(t1-t1[0]) + np.count_nonzero(t2-t2[0]) < 1:
            return (0<t1[0]<1) and (0<t2[0]<1)
    return False


def veclen(vec):
    return np.dot(vec,vec)**0.5


def nontriangulated_sides(segments, coordinates):
    # find out which of the edges contained in the list is not the side for 
    # two triangles. This happens if 
    # (1) the edge is on the boundary of the domain OR 
    # (2) there are edges missing in the domain, making an incomplete triangulation

    all_edge_segments = True
    hanging_segments = []
    segment_neighbor_count = [0 for s in segments]

    # create quick lookup of all edges emanating from each point
    node_nbrs = [[] for c in range(len(coordinates))]
    for si, s in enumerate(segments):
        node_nbrs[s[0]].append(si)
        node_nbrs[s[1]].append(si)

    # generate all possible triangles (a,b,c) where a,b,c are segments and the endpoints agree
    all_tris = set()
    for s0, s in enumerate(segments):
        for s1 in node_nbrs[s[0]]:
            if s1 != s0:
                for s2 in node_nbrs[s[1]]:
                    if (s2 != s1) and (s2 != s0):
                        if len(set([segments[si][j] for si in [s1,s2] for j in [0,1]])) == 3:
                            all_tris.add(tuple(sorted([s0,s1,s2])))
    # now allocate each triangle as a neighbor to its 3 edges
    for t in all_tris:
        for ti in t:
            segment_neighbor_count[ti] += 1

    # tabulate how many segments are missing at least one triangle neighbor
    for si, s in enumerate(segments):
        if segment_neighbor_count[si] < 2:
            hanging_segments.append(s)
            if not all([0 in coordinates[p] for p in s]):
                all_edge_segments = False

    return all_edge_segments, hanging_segments, segment_neighbor_count 


def all_cycles(segments):
    # find all cycles in an undirected graph defined by the list
    # of segments of the form [v1,v2] to denote an edge between vertices v1 and v2.
    # note: this routine ignores multi-edges, and it does not require the vertices 
    # to have index set {0,..., len(vertices)}.

    def dfs_all_cycles(segments, es_at, start_vertex, end_vertex):
        the_list = [(start_vertex, [])]
        while the_list:
            state, path = the_list.pop()
            if path and state == end_vertex:
                yield path
                continue
            for next_edge in es_at[state]:
                if next_edge in path:
                    continue
                next_state = [x for x in segments[next_edge] if x != state][0]
                the_list.append((next_state, path + [next_edge]))

    cycles = []
    added = []
    vertices = list(set([v for s in segments for v in s]))
    indexed_edges = [[vertices.index(v) for v in s] for s in segments]
    edges_out_of = [[] for i in range(len(vertices))]
    for ei, e in enumerate(indexed_edges):
        edges_out_of[e[0]].append(ei)
        edges_out_of[e[1]].append(ei)

    for node in range(len(vertices)):
        cycles_at_node = [path for path in dfs_all_cycles(indexed_edges, edges_out_of, node, node)]
        for c in cycles_at_node:
            if tuple(sorted(c)) not in added:
                added.append(tuple(sorted(c)))
                cycles.append(c)
    return cycles


def find_linear_groups(edge_indices, segments, coordinates):
    # This routine takes a set of edge_indices, and returns a list of
    # lists, where the edge indices are grouped into lists such that 
    # each list consists of edges that are aligned along a single line. 
    to_return = []
    edge_remaining = [True for e in range(len(segments))]

    for e in edge_indices:
        if edge_remaining[e]:
            current_group = [e]
            e_pts = [coordinates[s] for s in segments[e]]
            for other_e in edge_indices:
                if (other_e != e) and edge_remaining[other_e]:
                    if all([is_collinear(e_pts, coordinates[p]) for p in segments[other_e]]):
                        current_group.append(other_e)
                        edge_remaining[other_e] = False
            edge_remaining[e] = False
            to_return.append(current_group)
    return to_return


def identify_regular_triangles(segments, segment_adjacency_count, coordinates):
    # Take a set of segments and the number of neighboring triangles 
    # for each, and return a list of 'regular triangles' i.e. triangles
    # whose edges each consist of r segments, for some r>1.

    ts = []
    # make a list of the segments that already have 2 neighbors (i.e. they're done)
    completed_segments = [x > 1 for x in segment_adjacency_count]
    interior_segments = [iy for iy, y in enumerate(segments) if not (0 in coordinates[y[0]] and 0 in coordinates[y[1]])]

    # add in all the segments that are on an edge and already have a neighbor
    for si, s in enumerate(segments):
        if segment_adjacency_count[si] > 0 and (si not in interior_segments):
            completed_segments[si] = True

    remaining_segments = [tuple(s) for si, s in enumerate(segments) if not completed_segments[si]]
    remaining_segments_index = [si for si, s in enumerate(segments) if not completed_segments[si]]
    cycles = all_cycles(remaining_segments)

    if len(cycles) > 0:
        for cycle in cycles:
            indexed_cycle = [remaining_segments_index[e] for e in cycle]
            current_triangle = find_linear_groups(indexed_cycle, segments, coordinates)
            num_edges_per_side = [len(x) for x in current_triangle]

            # remove cycles that aren't triangle, cycles that aren't regular triangles, and the cycle that contains all the outer edges
            if all([x not in interior_segments for x in indexed_cycle]) \
               or len(current_triangle) != 3 \
               or (max(num_edges_per_side) - min(num_edges_per_side)) > 0:
                pass
            else:
                ts.append([[segments[s] for s in t] for t in current_triangle])
    return ts


def tesselate(triangle, coordinates):
    if not (len(triangle[0]) == len(triangle[1]) == len(triangle[2])):
        print("Error in tesselate routine: not a regular triangle.")
        return [[]],[]

    r = len(triangle[0]) - 1

    def order_segments(segs):
        [b,e] = segs[0]
        l = [tuple(segs[0])]
        met = [False for s in segs]
        met[0] = True

        while not all(met):
            for si, s in enumerate(segs):
                if not met[si]:
                    if e in s:
                        other = s[1-s.index(e)]
                        l.append((e, other))
                        e = other
                        met[si] = True
                    elif b in s:
                        other = s[1-s.index(b)]
                        l = [(other, b)] + l
                        b = other
                        met[si] = True
        return l

    # order the segments along each side so that side1 and side2 emanate from the 
    # same vertex, and side3 begins at the end of side1 and ends at the end of side2. 
    side1 = order_segments(triangle[0])
    side2 = order_segments(triangle[1])
    side3 = order_segments(triangle[2])
    if side3[0][0] == side1[0][0]:
        s2 = side2
        side2 = side3
        side3 = s2
    elif side3[-1][1] == side1[0][0]:
        s2 = side2
        side2 = side3
        side3 = s2
    if side2[-1][1] == side1[0][0]:
        side2 = [[y[1],y[0]] for y in side2[::-1]]
    if side3[-1][1] == side1[-1][1]:
        side3 = [[y[1],y[0]] for y in side3[::-1]]

    # generate new points:
    new_points = []
    for row in range(r - 1):
        num_pts_in_row = r - (row + 1)
        if num_pts_in_row > 1:
            row_endpoints = [coordinates[side2[row][1]], coordinates[side3[row][1]]]
            row_pts = list(zip(*[np.linspace(coordinates[0][i], coordinates[1][i], num_pts_in_row) for i in range(3)]))
            new_points.append(row_pts)
        elif num_pts_in_row > 0:
            row_endpoints = [coordinates[side2[row][1]], coordinates[side3[row][1]]]
            row_pts = [0.5*(row_endpoints[0][i]+row_endpoints[1][i]) for i in range(3)]
            new_points.append(row_pts)

    # create lookup table of index points for new_segments and new_points
    vertex_map = [s[0] for s in side1] + [side1[-1][1]]
    b = 0
    for row in range(r):
        vertex_map.append(side2[row][1])
        if row < (r - 1):
            vertex_map.extend([-(x+1) for x in range(b, b + r - (1 + row))])
            b += r - (1 + row)
        vertex_map.append(side3[row][1])
    vertex_map.append(side2[-1][1])

    # now create segments with appropriately indexed new_pts values
    new_segments = []
    offset = 0
    for row in range(r):
        new_segments.extend([[vertex_map[offset+c], vertex_map[offset+(r+2-row)+c+j]] for j in [-1,0] for c in range(1, r+1-row)])
        new_segments.extend([[vertex_map[offset+(r+2-row)+c+j] for j in [-1,0]] for c in range(1, r+1-row)])
        offset += r + 2 - row

    return new_points, new_segments


def generate_initial_rays(R, eis, Li):
    all_rays = []
    strengths = []
    for i in range(3):

        ei = tuple([0 if x != i else R for x in range(3)])
        deleted_lattice_points = [tuple(eis[(i+1)%3]), tuple(eis[(i+2)%3])] + [x for x in Li if R not in tuple(x)]

        indices_of_cvxhull_pts = convex_hull(deleted_lattice_points)
        # extract the coordinatens of the points on the convex hull
        pts = [deleted_lattice_points[k] for k in indices_of_cvxhull_pts] 

        # convex_hull includes all lattice points along sides of the simplex, so we 
        # need to remove all the edge-lattice points except the two closest ones that lie along the 
        # sides emanating from ei 
        side_pt_indices = [[j for j,p in enumerate(pts) if p[(i+k)%3] == 0] for k in range(3)]
        side_1_closest_pt = pts[sorted([(veclen(np.array(pts[pi]) - np.array(ei)), pi) for pi in side_pt_indices[1]])[0][1]]
        side_2_closest_pt = pts[sorted([(veclen(np.array(pts[pi]) - np.array(ei)), pi) for pi in side_pt_indices[2]])[0][1]]
        side_pts = list(set([y for x in side_pt_indices for y in x]))
        pts = [side_1_closest_pt] + [x for i, x in enumerate(pts) if i not in (side_pts)] + [side_2_closest_pt]

        Lis = [x - eis[i] for x in pts]

        # generate initial rays and strengths for e_i
        for j, x in enumerate(ordered_rays(Lis)):
            if x[1] > 0:
                all_rays.append([ei, tuple(pts[x[0]])])
                strengths.append(int(x[1]))

    # generate index-valued copies of rays and points so referencing can be done without coordinates
    return all_rays, strengths


def Reids_recipe(segments, strengths, not_done, potential_segments, longest_extension, coordinates):

    new_segments = [(-1, s[0], s[1], strengths[si]) for si, s in enumerate(segments)]
    children = [[] for si in range(len(segments))]

    # make sure each initial segment that meets with an initial segment of greater strength 
    # is marked as "finished" so that it is not extended further into the domain. 
    not_finished = [True for s in range(len(strengths))]
    for i, si in enumerate(segments):
        segs_at = [j for j, sj in enumerate(segments) if sj[1] == si[1]]
        strengths_at = [strengths[j] for j in segs_at]
        for j in segs_at:
            not_finished[j] = False
            if (strengths[j] == max(strengths_at)) and (strengths_at.count(max(strengths_at)) < 2):
                not_finished[j] = True

    # add in extensions to each remaining non-finished segment
    for si, s in enumerate(segments):
        if not_finished[si]:
            extensions_for_s = [(ps[1], ps[2], 0) for ps in potential_segments if ps[0] == si]

            if len(extensions_for_s) > 0:
                indices_from = [si] + [len(new_segments) + j for j in range(len(extensions_for_s) - 1)]
                extensions_for_s = [(indices_from[i],e[0],e[1],e[2]) for i, e in enumerate(extensions_for_s)]

                added_indices = [len(new_segments) + j for j in range(len(extensions_for_s))]
                children[si] = added_indices
                for e in range(len(extensions_for_s) - 1):
                    children.append(added_indices[e+1:])
                children.append([])
                new_segments.extend(extensions_for_s)

    segments = new_segments

    segs_at_pt = [[] for c in coordinates]
    for si, s in enumerate(segments):
        segs_at_pt[s[2]].append(si)

    still_updating = True
    ctr = 0
    while (still_updating and (ctr < longest_extension*len(segments))):
        ctr += 1
        still_updating = False

        # order segments by their strengths
        current_segs = sorted([(segments[si][-1], si) for si in range(len(segments))])

        for (seg_i_strength, seg_i) in current_segs[::-1]:
            segment = segments[seg_i]

            # check if this segment has nonzero strength, and that it has possible extensions
            if segment[-1] >= 0 and len(children[seg_i]) > 0:

                segments_at_endpoint = [s for s in segs_at_pt[segment[2]] if segments[s][-1] >= 0]
                strengths_at_endpoint = [segments[s][-1] for s in segments_at_endpoint]

                first_child = children[seg_i][0]

                if segment[-1] < max(strengths_at_endpoint):
                    if segments[first_child][-1] >= 0:
                        still_updating = True

                    for c in children[seg_i]:
                        segments[c] = (segments[c][0], segments[c][1], segments[c][2], -1)
                else:
                    if strengths_at_endpoint.count(max(strengths_at_endpoint)) < 2:

                        not_intersects = True
                        for sj, s in enumerate(segments):
                            if sj != first_child and s[-1] > 0:
                                if intersects([segments[sj][1], segments[sj][2]], [segments[first_child][1], segments[first_child][2]], coordinates):
                                    not_intersects = False
                                    break
                        if not_intersects:
                            child_strength = segment[-1] + 1 - len(segments_at_endpoint)
                            if child_strength != segments[first_child][-1]:
                                 segments[first_child] = (segments[first_child][0], segments[first_child][1], segments[first_child][2], child_strength)
                                 still_updating = True
                        else:
                            if segments[first_child][-1] >= 0:
                                still_updating = True
                            for c in children[seg_i]:
                                segments[c] = (segments[c][0], segments[c][1], segments[c][2], -1)

            if still_updating:
                break

    strengths = [s[-1] for s in segments if s[-1] > 0]
    segments = [[s[1],s[2]] for si, s in enumerate(segments) if s[-1] > 0]

    return segments


def triangulation(R,a,b,c):
    # vertices defining the junior simplex (scaled by length R so we deal with integers instead of fractions)
    eis = [np.array([0 if j != i else R for j in range(3)]) for i in range(3)]

    # interior lattice points for junior simplex (also scaled by r)
    coordinates = lattice_points_nonunit(R,a,b,c) # list of arrays
    tuplePts = [tuple(x) for x in coordinates] # list of tuples (used for finding index of point)

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #       STEP 1: GENERATE INITIAL RAYS AND THEIR STRENGTHS       #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # create rays L_{i,j} emanating from each corner e_i, each with 
    # its associated 'strength' (from the Jung-Hirzebruch relation)
    segments_with_coords, strengths = generate_initial_rays(R, eis, coordinates)
    segments = [tuple([tuplePts.index(r[j]) for j in range(2)]) for r in segments_with_coords]


    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #       STEP 2: CALCULATE POSSIBLE EXTENSIONS FOR EACH RAY      #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # not_done and potential_segments keep track of which rays can be "extended" past their endpoints (and how far)
    not_done = [True if strengths[si] != 0 else False for si, s in enumerate(segments)]
    potential_segments = []
    longest_extension = 0
    for ir, r in enumerate(segments_with_coords):
        if strengths[ir] > 0:
            pts_along_r = [(0, segments[ir][1])] + sorted([(veclen(p-np.array(r[1])), pi) for pi, p in enumerate(coordinates) if ((pi not in segments[ir]) and (is_collinear(r, p)))])
            potential_segments.extend([[ir, pts_along_r[j][1], pts_along_r[j+1][1]] for j in range(len(pts_along_r) - 1)])
            longest_extension = max(longest_extension, (len(pts_along_r)))


    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #       STEP 3: FIND INTERSECTIONS OF RAYS/EXTENSIONS           #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # now keep track of main list of segments that are going to be added to create a triangulation.
    segments = Reids_recipe(segments, strengths, not_done, potential_segments, longest_extension, coordinates)
    # add segments that lie along the sides
    for i in range(3):
        points_along_side = sorted([tuple(x) for x in coordinates if x[i] == 0])
        segments_along_side = [[tuplePts.index(points_along_side[ip]), tuplePts.index(points_along_side[ip+1])] for ip in range(len(points_along_side)-1)]
        segments.extend(segments_along_side)
        strengths.extend([0 for x in range(len(segments_along_side))])
    segments = list(set([tuple(sorted(s)) for s in segments]))


    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #     STEP 4: TESSELATE REMAINING r-SIDED TRIANGLES INTO r^2    #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    all_edges, hanging_segments, segment_neighbor_count = nontriangulated_sides(segments, coordinates)
    if not all_edges:
        regular_triangles = identify_regular_triangles(segments, segment_neighbor_count, coordinates)

        for t in regular_triangles:
            extra_points, extra_segments = tesselate(t, coordinates)

            for e_s in extra_segments:
                reindexed_segment = tuple([e_s[i] if e_s[i] >= 0 else (len(coordinates) - (e_s[i] + 1)) for i in range(2)])
                segments.append(reindexed_segment)

            if len(extra_points) > 0:
                coordinates = coordinates + extra_points

    return segments, coordinates


def rotateSimplexToPlane(coords):
    n = np.array([1,1,1])
    n = n / np.sqrt(3)
    m1 = np.matrix([
        [ n[0] / (n[0]**2 + n[1]**2), n[1] / (n[0]**2 + n[1]**2), 0], 
        [-n[1] / (n[0]**2 + n[1]**2), n[0] / (n[0]**2 + n[1]**2), 0], 
        [0, 0, 1]])
    n = np.matmul(n, m1.transpose()).tolist()[0]
    m2 = np.matrix([
        [n[2],0,-n[0]],
        [0,1,0],
        [n[0],0, n[2]]])
    return [np.matmul(np.matrix(c), np.matmul(m2,m1).transpose()).tolist()[0] for c in coords]


def curve_type(edges, triangles, edge_to_triangle, coordinates, e):
    if isinstance(e, list):
        if isinstance(edges[0], list):
            edges = [tuple(i) for i in edges]
        ei = edges.index(tuple(e))
    else:
        ei, e = e, edges[e]

    e2t = sorted(edge_to_triangle[ei])
    if len(e2t) < 2: # check if its an edge on the boundary
        return (0,0)

    # unpack the vertices on the edge and those it would be flopped to 
    vertices = set([tuple(coordinates[i]) for t in e2t for s in triangles[t] for i in edges[s]])
    # vertices of existing edge
    e_verts = [tuple(coordinates[e[0]]), tuple(coordinates[e[1]])]

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
    if abs(angle1) < tol or abs(angle2) < tol:
        return (-2, 0)
    if np.sign(angle1) == np.sign(angle2):
        # check that neither angle is zero
        if (abs(angle1) > tol) and (abs(angle2) > tol):
            return (-1,-1)
        return (-2,0)
    return (-3,0)

def QPFromTriangulation(R,a,b,c):
    (edges, coordinates) = triangulation(R,a,b,c)
    extremal_edges = [ei for ei, e in enumerate(edges) if (R in tuple(coordinates[e[0]]) and R in tuple(coordinates[e[1]]))]
    coordinates = [[c[0], c[1]] for c in rotateSimplexToPlane(coordinates)]

    tri = make_triangulation(edges, extremal_edges)
    triangles = tri[0]
    edge_to_triangle = tri[1]

    QP_edges = [(t[i%3], t[(i+j)%3]) \
            for j in [-1,1] for i in range(1, 4) \
            for t in triangles \
            if ((len(edge_to_triangle[t[i%3]]) > 1) and (len(edge_to_triangle[t[(i+j)%3]]) > 1))]

    for i, e in enumerate(edges):
        if len(edge_to_triangle[i]) > 1:
            if curve_type(edges, triangles, edge_to_triangle, coordinates, list(e)) == (-2,0):
                QP_edges.append((i,i))
            if curve_type(edges, triangles, edge_to_triangle, coordinates, list(e)) == (-3,0):
                QP_edges.append((i,i))
                QP_edges.append((i,i))

    e_reorder = [i for i in range(len(edges)) if len(edge_to_triangle[i]) > 1]
    QP_edges = [(e_reorder.index(x[0]), e_reorder.index(x[1])) for x in QP_edges]
    p = {}

    positions = [[0.5*(coordinates[e[0]][0]
                      +coordinates[e[1]][0]), 
                  0.5*(coordinates[e[0]][1]
                      +coordinates[e[1]][1])] for ei, e in enumerate(edges) if len(edge_to_triangle[ei]) > 1]
    return QuiverWithPotential(QP_edges, positions=positions)

#R,a,b,c=11,1,2,8
#R,a,b,c=30,25,2,3
#R,a,b,c=6,1,2,3
R,a,b,c = 25,1,3,21

QP = QPFromTriangulation(R,a,b,c)
QP.toJSON("%d_%d_%d_%d.JSON"%(R,a,b,c))

