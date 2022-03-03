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
    numerator = np.array(pt) - np.array(ray[0])
    denominator = np.array(ray[1]) - np.array(ray[0])

    if set(numerator.nonzero()[0]) != set(denominator.nonzero()[0]):
        return False

    frac = numerator / denominator
    return frac.max() - frac.min() < tol


def intersects(seg1, seg2, pts):
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
                            not_intersects = all([~intersects([ps[1], ps[2]], s, points) for s in all_segments])
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

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for s in segments:
            xp, yp, zp = np.matrix([points[s[0]], points[s[1]]]).transpose().tolist()
            ax.plot3D(xp,yp,zp)
        plt.show()

    # add segments that lie along the sides
    for i in range(3):
        points_along_side = sorted([tuple(x) for x in points if x[i] == 0])
        segments_along_side = [[tuplePts.index(points_along_side[ip]), tuplePts.index(points_along_side[ip+1])] for ip in range(len(points_along_side)-1)]
        segments.extend(segments_along_side)
        strengths.extend([0 for x in range(len(segments_along_side))])
        

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    #     STEP 4: TESSELATE REMAINING r-SIDED TRIANGLES INTO r^2    #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for s in segments:
        xp, yp, zp = np.matrix([points[s[0]], points[s[1]]]).transpose().tolist()
        ax.plot3D(xp,yp,zp)
    plt.show()



R,a,b,c=6,1,2,3
R,a,b,c=30,25,2,3
#R,a,b,c=11,1,2,8
triangulation(R,a,b,c)
