import numpy as np
import copy
import matplotlib.pyplot as plt
from mutations import *


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
            row_pts = list(zip(*[np.linspace(row_endpoints[0][i], row_endpoints[1][i], num_pts_in_row+2) for i in range(3)]))
            new_points.extend([tuple(x) for x in row_pts[1:-1]])
        elif num_pts_in_row > 0:
            row_endpoints = [coordinates[side2[row][1]], coordinates[side3[row][1]]]
            row_pts = tuple([0.5*(row_endpoints[0][i]+row_endpoints[1][i]) for i in range(3)])
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

    return new_points, list(set([tuple([min(x), max(x)]) for x in new_segments]))


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


def Reids_recipe(segments, strengths, potential_segments, longest_extension, coordinates):

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
    # potential_segments keeps track of which rays can be "extended" past their endpoints (and how far)
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
    segments = Reids_recipe(segments, strengths, potential_segments, longest_extension, coordinates)
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
    if len(o_verts) > 1:
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
    return (0,0)


def drawTriangulation(t):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for s in t[0]:
        xp, yp, zp = np.matrix([t[1][s[0]], t[1][s[1]]]).transpose().tolist()
        ax.plot3D(xp,yp,zp)
    ax.view_init(elev=45., azim=45)
    plt.show()


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
R,a,b,c = 60,49,6,5

t = triangulation(R,a,b,c)
drawTriangulation(t)
QP = QPFromTriangulation(R,a,b,c)
QP.toJSON("%d_%d_%d_%d.JSON"%(R,a,b,c))

