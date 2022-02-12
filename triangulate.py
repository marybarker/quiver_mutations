import numpy as np
import copy
import random
from itertools import permutations
import matplotlib.pyplot as plt


def interior_lattice_points_nonunit(r=6,a=1,b=2,c=3):
    # Lattice points in the junior simplex lie in the plane x+y+z=1 and are of the form 
    # (i/r,j/r,k/r) where i+j+k=r. 
    # This routine returns the triples (i,j,k) instead of scaling by 1/r. 
    points = []
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
    l0,lk = L[-1],L[0]
    iLs = L[1:-1]
    if len(iLs) < 2:
        return [(0, 1), (1, lk[0] + l0[0]), (2, 1)]
    else:
        # find L1 and L2 as a starting point
        firstPair = []
        for iL1,L1 in enumerate(iLs):
            if len(firstPair) < 1:
                for iL2,L2 in enumerate(iLs):
                    if iL2 != iL1:
                        vals = l0 + L2
                        if set(vals.nonzero()[0]) == set(L1.nonzero()[0]):
                            vals = vals / L1
                            if (vals.max() - vals.min()) < 1.0e-6:
                                firstPair = [iL1, iL2, vals[vals.nonzero()[0][0]]]
                                break
        if len(firstPair) < 1:
            print("ERROR IN ORDERED_RAYS ROUTINE!!!\nNever found a first pair")
            exit(1)
        toRet = [(0, 1), (firstPair[0], firstPair[2])]

        met = [firstPair[0]]
        L0 = iLs[firstPair[0]]
        iL1 = firstPair[1]
        L1 = iLs[iL1]

        for i in range(len(iLs) - 1):
            not_found = True
            for iL2, L2 in enumerate(iLs):
                if not_found:
                    if (iL2 not in met) and iL1 != iL2:
                        vals = L0 + L2
                        if set(vals.nonzero()[0]) == set(L1.nonzero()[0]):
                            vals = vals / L1

                            if (vals.max() - vals.min()) < 1.0e-6:
                                not_found = False
                                toRet.append((iL1, vals[vals.nonzero()[0][0]]))
                                met.append(iL1)

                                # now update pointer to next step
                                L0 = iLs[iL1]
                                iL1 = iL2
                                L1 = iLs[iL1]

        nz = L1.nonzero()[0][0]
        # add second-to-last ray (its L2 will be the last ray)
        toRet.append((iL1, (L0[nz] + lk[nz]) / L1[nz]))
        # add the last ray with a strength of 1
        toRet.append([len(L)-1, 1])
        return toRet


def is_collinear(ray, pt, tol=1.0e-6):
    numerator = np.array(pt) - np.array(ray[0])
    denominator = np.array(ray[1]) - np.array(ray[0])

    if set(numerator.nonzero()[0]) != set(denominator.nonzero()[0]):
        return False

    frac = numerator / denominator
    return frac.max() - frac.min() < tol


def triangulation(R,a,b,c):
    # vertices defining the junior simplex (scaled by length R so we deal with integers instead of fractions)
    eis = [np.array([0 if j != i else R for j in range(3)]) for i in range(3)]

    # interior lattice points for junior simplex
    Li = interior_lattice_points_nonunit(R,a,b,c)

    # create rays L_{i,j} emanating from each corner e_i, each with 
    # its associated 'strength' (from the Jung-Hirzebruch relation)
    all_rays = []
    strengths = []
    for i in range(3):
        # get the points of the convex hull of {Delta + interior_points}\{e_i}
        points = [eis[(i+2)%3], eis[(i+1)%3]] + Li
        indices_of_cvxhull_pts = convex_hull(points)
        pts = [points[k] for k in indices_of_cvxhull_pts]
        # reindex so that first and last points of Li correspond to the edges of the sinmplex
        pts = [pts[0]] + pts[2:] + [pts[1]]

        side1_pts = sorted([((x-eis[i]).sum(), ix) for ix, x in enumerate(points) if is_collinear(eis[(i+2)%3]-eis[i], x)])
        if len(side1_pts) > 0:
            side1_pts = [y[1] for y in side1_pts]
            pts = [pts[side1_pts[0]]] + [x for ix, x in enumerate(points) if ix not in side1_pts] 

        side2_pts = sorted([((x-eis[i]).sum(), ix) for ix, x in enumerate(points) if is_collinear(eis[(i+1)%3]-eis[i], x)])
        if len(side2_pts) > 0:
            side2_pts = [y[1] for y in side2_pts]
            pts = [x for ix, x in enumerate(points) if ix not in side2_pts] + [pts[side2_pts[0]]]
    
        # now make these into rays emanating from e_i
        Lis = [x - eis[i] for x in pts]
    
        # generate initial rays and strengths for e_i
        for j, x in enumerate(ordered_rays(Lis)):
            all_rays.append([eis[i], np.array(pts[j])])
            strengths.append(int(x[1]))

    # generate index-valued copies of rays and points so referencing can be done without coordinates
    points = eis + Li
    tuplePts = [tuple(x) for x in points]
    segments = [[tuplePts.index(tuple(r[j])) for j in range(2)] for r in all_rays]

    # not_done and potential_segments keep track of which rays can be "extended" past their endpoints (and how far)
    not_done = [True for s in segments]
    potential_segments = []
    longest_extension = 0
    for ir, r in enumerate(all_rays):
        if strengths[ir] > 1:
            pts_along_r = [(0, segments[ir][1])] + sorted([((p-r[1]).sum(), pi) for pi, p in enumerate(points) if ((pi not in segments[ir]) and (is_collinear(r, p)))])
            potential_segments.extend([[ir, pts_along_r[i][1], pts_along_r[i+1][1]] for i in range(len(pts_along_r) - 1)])
            longest_extension = max(longest_extension, (len(pts_along_r)))

    # now keep track of main list of segments that are going to be added to create a triangulation.
    all_segments = copy.deepcopy(segments)
    all_strengths = copy.deepcopy(strengths)
    for extensions in range(longest_extension + 1):
        for i in range(len(segments)):

            # if the endpoint of this segment has not had an intersection computed already
            if not_done[i] and (strengths[i] > 1):
                pt = segments[i][1] # extract the endpoint and all segments that hit that point
                segments_at_pt = [j for j,s in enumerate(segments) if (s[1] == pt and strengths[j] > 1 and not_done[j])]

                # if we have multiple segments intersecting at the point:
                if len(segments_at_pt) > 1:

                    for j in segments_at_pt:
                        not_done[j] = False

                    strengths_at_pt = [strengths[j] for j in segments_at_pt]
                    mx = max(strengths_at_pt)
                    # first check if all rays meeting at the point have equal strength (or there are multiple 'winners')
                    if strengths_at_pt.count(mx) < 2:
                        ij = strengths_at_pt.index(mx)
                        j = segments_at_pt[ij]

                        new_s = mx - (len(strengths_at_pt) - 1)
                        ps_from_j = [psi for psi, ps in enumerate(potential_segments) if (ps[0] == j)]
                        ps = potential_segments[ps_from_j[0]]

                        all_segments.append([ps[1], ps[2]])
                        all_strengths.append(new_s)
                        not_done.append(True)

                        # now update potential segments so that they emanate from the newly added segment
                        for k in ps_from_j[1:]:
                            oldps = potential_segments[k]
                            potential_segments[k] = [len(all_strengths) - 1, oldps[1], oldps[2]]

        segments = copy.deepcopy(all_segments)
        strengths = copy.deepcopy(all_strengths)


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for s in all_segments:
        xp, yp, zp = np.matrix([points[s[0]], points[s[1]]]).transpose().tolist()
        ax.plot3D(xp,yp,zp)
    plt.show()
    plt.clf()



R,a,b,c=13,1,5,7#11,1,2,8
triangulation(R,a,b,c)
