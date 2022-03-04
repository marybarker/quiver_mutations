import numpy as np
import copy
import random
from itertools import permutations
import matplotlib.pyplot as plt


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
        toRet = [(0,0)] + [(angles[i-1]+2, (oLs[i-1][np.flatnonzero(oLs[i])[0]] + oLs[i+1][np.flatnonzero(oLs[i])[0]]) / oLs[i][np.flatnonzero(oLs[i])[0]]) for i in range(1, len(oLs)-1)] + [(1,0)]
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
                    if all([is_collinear(e_pts, coordinates[p]) for p in segments[other_e] if p not in segments[e]]):
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
    elif side3[-1][0] == side1[0][0]:
        s3 = side2
        side2 = [[y[1],y[0]] for y in side3[::-1]]
        side3 = s3
    elif side2[-1][0] == side1[0][0]:
        side2 = [[y[1],y[0]] for y in side2[::-1]]
    if side3[-1][0] == side1[0][0]:
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
            row_pts = [0.5*(coordinates[0][i]+coordinates[1][i]) for i in range(3)]
            new_points.append(row_pts)

    # create lookup table of index points for new_segments and new_points
    vertex_map = [s[0] for s in side1] + [side1[-1][1]]
    b = 0
    for row in range(r):
        vertex_map.append(side2[row][1])
        vertex_map.extend([-(x+1) for x in range(b, b + r - 1 - row)])
        vertex_map.append(side3[row][1])
        b += r - 1 - row
    vertex_map.append(side2[-1][1])

    # now create segments with appropriately indexed new_pts values
    new_segments = []
    offset = 0
    for row in range(r):
        new_segments.extend([[vertex_map[offset+c], vertex_map[offset+(r+2-row)+c+j]] for j in [-1,0] for c in range(1, r+1-row)])
        new_segments.extend([[vertex_map[offset+(r+2-row)+c+j] for j in [-1,0]] for c in range(1, r+1-row)])
        offset += r + 2 - row

    return new_points, new_segments


def triangulation(R,a,b,c):
    # vertices defining the junior simplex (scaled by length R so we deal with integers instead of fractions)
    eis = [np.array([0 if j != i else R for j in range(3)]) for i in range(3)]

    # interior lattice points for junior simplex (also scaled by r)
    Li = lattice_points_nonunit(R,a,b,c)

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #       STEP 1: GENERATE INITIAL RAYS AND THEIR STRENGTHS       #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # create rays L_{i,j} emanating from each corner e_i, each with 
    # its associated 'strength' (from the Jung-Hirzebruch relation)
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
                all_rays.append([ei, tuple(pts[j])])
                strengths.append(int(x[1]))

    # generate index-valued copies of rays and points so referencing can be done without coordinates
    points = Li
    tuplePts = [tuple(x) for x in points]
    segments = [tuple([tuplePts.index(r[j]) for j in range(2)]) for r in all_rays]

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #       STEP 2: CALCULATE POSSIBLE EXTENSIONS FOR EACH RAY      #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # not_done and potential_segments keep track of which rays can be "extended" past their endpoints (and how far)
    not_done = [True if strengths[si] != 0 else False for si, s in enumerate(segments)]
    potential_segments = []
    longest_extension = 0
    for ir, r in enumerate(all_rays):
        if strengths[ir] > 0:
            pts_along_r = [(0, segments[ir][1])] + sorted([(veclen(p-np.array(r[1])), pi) for pi, p in enumerate(points) if ((pi not in segments[ir]) and (is_collinear(r, p)))])
            potential_segments.extend([[ir, pts_along_r[j][1], pts_along_r[j+1][1]] for j in range(len(pts_along_r) - 1)])
            longest_extension = max(longest_extension, (len(pts_along_r)))

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #       STEP 3: FIND INTERSECTIONS OF RAYS/EXTENSIONS           #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # now keep track of main list of segments that are going to be added to create a triangulation.
    all_segments = copy.deepcopy(segments)
    all_strengths = copy.deepcopy(strengths)
    for extensions in range(longest_extension):

        added_segments = False
        for i in range(len(segments)):

            # if the endpoint of this segment has not had an intersection computed already
            if not_done[i] and (strengths[i] > 0):
                pt = segments[i][1] # extract the endpoint and all segments that hit that point
                segments_at_pt = [j for j,s in enumerate(all_segments) if (s[1] == pt and all_strengths[j] > 0)]

                # if we have multiple segments intersecting at the point:
                if len(segments_at_pt) > 1:

                    for j in segments_at_pt:
                        not_done[j] = False

                    strengths_at_pt = [all_strengths[j] for j in segments_at_pt]
                    mx = max(strengths_at_pt)
                    # first check if all rays meeting at the point have equal strength (or there are multiple 'winners')
                    if strengths_at_pt.count(mx) < 2:
                        ij = strengths_at_pt.index(mx)
                        j = segments_at_pt[ij]

                        new_s = mx - (len(strengths_at_pt) - 1)
                        ps_from_j = [psi for psi, ps in enumerate(potential_segments) if (ps[0] == j)]
                        if len(ps_from_j) > 0:
                            ps = potential_segments[ps_from_j[0]]

                            all_segments.append((ps[1], ps[2]))
                            all_strengths.append(new_s)
                            not_done.append(True)
                            added_segments = True

                            # now update potential segments so that they emanate from the newly added segment
                            for k in ps_from_j[1:]:
                                oldps = potential_segments[k]
                                potential_segments[k] = [len(all_strengths) - 1, oldps[1], oldps[2]]
                            potential_segments = potential_segments[:ps_from_j[0]] + potential_segments[ps_from_j[0]+1:]
        if added_segments:
            segments = copy.deepcopy(all_segments)
            strengths = copy.deepcopy(all_strengths)

        # now double-check if there were no intersections (all segments need to be extended in order to reach another)
        if (not added_segments):
            for i in range(len(segments)):
                if strengths[i] > 0 and not_done[i]:
                    # get endpoint of segment
                    pt = segments[i][1]

                    # get the next 'potential' segment that extends the current one at the endpoint
                    ps_from_pt = [psi for psi, ps in enumerate(potential_segments) if (ps[0] == i)]

                    if len(ps_from_pt) > 0:
                        ps = potential_segments[ps_from_pt[0]]

                        # make sure that we only extend the current segment till it intersects with another
                        # i.e. only extend this segment if it doesn't intersect with another segment at its endpoint
                        other_potential_segments_hitting_pt = len([1 for x in potential_segments if x[2] == ps[2]])
                        previous_segments_hitting_pt = len([1 for ix, x in enumerate(segments) if x[1] == pt and i != ix])

                        if (other_potential_segments_hitting_pt + previous_segments_hitting_pt) < 3:

                            # double-check if this potential segment intersects a previously existing one
                            not_intersects = all([not intersects([ps[1], ps[2]], s, points) for s in all_segments])
                            if not_done[ps[0]] and not_intersects:
                                all_segments.append((ps[1], ps[2]))
                                all_strengths.append(all_strengths[ps[0]])
                                not_done.append(True)
                                potential_segments[ps_from_pt[0]] = [len(all_strengths) - 1, ps[1], ps[2]]

                            for psi, ps in enumerate(potential_segments):
                                if ps[0] == i:
                                    potential_segments[psi] = [len(all_strengths) - 1, ps[1], ps[2]]
                            potential_segments = potential_segments[:ps_from_pt[0]] + potential_segments[ps_from_pt[0]+1:]
                            not_done[i] = False
            segments = copy.deepcopy(all_segments)
            strengths = copy.deepcopy(all_strengths)

        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #for s in segments:
        #    xp, yp, zp = np.matrix([points[s[0]], points[s[1]]]).transpose().tolist()
        #    ax.plot3D(xp,yp,zp)
        #plt.show()

    # add segments that lie along the sides
    for i in range(3):
        points_along_side = sorted([tuple(x) for x in points if x[i] == 0])
        segments_along_side = [[tuplePts.index(points_along_side[ip]), tuplePts.index(points_along_side[ip+1])] for ip in range(len(points_along_side)-1)]
        segments.extend(segments_along_side)
        strengths.extend([0 for x in range(len(segments_along_side))])
    segments = list(set([tuple(sorted(s)) for s in segments]))

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #     STEP 4: TESSELATE REMAINING r-SIDED TRIANGLES INTO r^2    #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #for si, s in enumerate(segments):
    #    print("segment %d : "%si + str(points[s[0]]) + ", "+str(points[s[1]]))
    all_edges, hanging_segments, segment_neighbor_count = nontriangulated_sides(segments, points)
    for s in sorted(segments):
        print(s)
    if not all_edges:
        regular_triangles = identify_regular_triangles(segments, segment_neighbor_count, points)

        print("There are %d regular triangles"%len(regular_triangles))
        for t in regular_triangles:
            extra_points, extra_segments = tesselate(t, points)
            print("adding in %d extra points and %d extra segments"%(len(extra_points), len(extra_segments)))

            for e_s in extra_segments:
                reindexed_segment = tuple([e_s[i] if e_s[i] >= 0 else (len(points) - (e_s[i] + 1)) for i in range(2)])
                segments.append(reindexed_segment)

            if len(extra_points) > 0:
                points = points + extra_points

    print(points, segments)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for s in segments:
        xp, yp, zp = np.matrix([points[s[0]], points[s[1]]]).transpose().tolist()
        ax.plot3D(xp,yp,zp)
    plt.show()



R,a,b,c=6,1,2,3
#R,a,b,c=30,25,2,3
#R,a,b,c=11,1,2,8
triangulation(R,a,b,c)
