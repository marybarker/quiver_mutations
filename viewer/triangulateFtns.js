var globalTriangulation = [
	[[0,5],[5,4],[4,1],[1,2],[2,6],[6,0],
	[6,5],[5,3],[4,3],[1,3],[2,3],[6,3]], 
	[[600,0,0],[0,600,0],[0,0,600],
	[100,200,300],[200,400,0],[400,200,0],[300,0,300]]
];

function allCycles(segments) {
    // find all cycles in an undirected graph defined by the list of segments 
    // of the form [v1,v2] to denote an edge between the vertices v1 and v2.
    // note: this routine ignores multi-edges, and it does not require the vertices 
    // to have index set {0,..., vertices.length}
    var cycles = [];
    var added = [];
    var vertices = new Set(segments.reduce(function(a, b) {return a.concat(b);}));
    vertices = Array.from(vertices);
    var indexed_edges = segments.map(s => s.map(x => vertices.indexOf(x)));
    var edges_out_of = vertices.map(x => []);

    for (let ei = 0; ei < indexed_edges.length; ei++) {
        let e = indexed_edges[ei];
        edges_out_of[e[0]].push(ei);
        edges_out_of[e[1]].push(ei);
    }

    for (let node = 0; node < vertices.length; node++) {
        let cycles_at_node = dfsAllCycles(indexed_edges, edges_out_of, node, node);
        for (let ci = 0; ci < cycles_at_node.length; ci++) {
            let c = cycles_at_node[ci];
            sortedc = [...c].sort().toString();
            if (!added.includes(sortedc)) {
                cycles.push(c);
                added.push(sortedc);
            }
        }
    }
    return cycles;
}

function argsort(arr) {
    // dsu algorithm for argsort
    let decor = (v, i) => [v, i];
    let undecor = a => a[1];
    return arr.map(decor).sort().map(undecor);
}

function callTriangulation() {
    var a = parseFloat(document.getElementById("A_value").value);
    var b = parseFloat(document.getElementById("B_value").value);
    var c = parseFloat(document.getElementById("C_value").value);
    var r = a + b + c;
    var t = triangulation(r,a,b,c);
    globalTriangulation = t;
    viewTriangulation();
}

function convexHull(planar_points) {
    // performs a quickhull search to find which points in the set of 
    // planar_points form the convex hull.
    // This routine assumes that the set of planar_points are points 
    // in 3D, but which lie in some plane. Also assumes that the first 
    // two points in the set of planar_points forms an edge of the 
    // convex hull. The more general version defines the initial starting 
    // points p0 and p1 as any two points with maximal distance between them
    var p0 = planar_points[0];
    var p1 = planar_points[1];
    var s1 = findHull(planar_points.slice(2), p0, p1);
    var s2 = findHull(planar_points.slice(2), p1, p0);
    var points = unique(s1.concat(s2)).map(p => planar_points.indexOf(p));
    return [...points].sort().map(p => planar_points[p]);
}

function count(l, v) {
    return l.filter(x => x == v).length;
}

function crossProduct(v1, v2) {
    return [v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]];
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
function curveType(edges, triangles, edge_to_triangle, coordinates, e) {
    return [0,0];
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

function dfsAllCycles(segments, es_at, start_vertex, end_vertex) {
    // helper function for allCycles routine
    var toRet = [];
    var the_list = [[start_vertex, []]];
    while (the_list.length > 0) {
        var output = the_list.pop();
        var state = output[0];
        var path = output[1];
        if ((path.length > 0) && (state == end_vertex)) {
            toRet.push(path);
            continue;
        } 
        for (let nei = 0; nei < es_at[state].length; nei++) {
            let next_edge = es_at[state][nei];
            if (path.includes(next_edge)) {
                continue;
            }
            var next_state = segments[next_edge].filter(x => x != state)[0];
            the_list.push([next_state, path.concat(next_edge)]);
        }
    }
    return toRet;
}

function dotProduct(v1, v2) {
    return v1.map((x, i) => v1[i] * v2[i]).reduce((m, n) => m + n);
}

function findHull(point_set, a, b) {
    // line segment from a to b and its length
    var AB = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    var rays = point_set.map(p => [p[0] - a[0], p[1] - a[1], p[2] - a[2]]);

    // find all the points on one side of the segment ab
    var pts = [];
    for (let i = 0; i < rays.length; i++) {
        if (crossProduct(AB, rays[i]).filter(x => x < 0).length < 1) {
            pts.push(point_set[i]);
        }
    }
    rays = pts.map(p => [p[0] - a[0], p[1] - a[1], p[2] - a[2]]);

    if (pts.length > 0) {
        // calculate the length of semgent ab
        var ABL = dotProduct(AB, AB);
        // lenght of vector from each pt to the line segment ab
        var distances = rays.map(r => [
		r[0] - (dotProduct(AB, r)/ABL)*AB[0],
		r[1] - (dotProduct(AB, r)/ABL)*AB[1],
		r[2] - (dotProduct(AB, r)/ABL)*AB[2]]
            ).map(x => dotProduct(x, x));
        var maxPt = distances.indexOf(Math.max(...distances));
        var dpts = pts.slice(0,maxPt).concat(pts.slice(maxPt+1));
        var s1 = findHull(dpts, a, pts[maxPt]);
        var s2 = findHull(dpts, pts[maxPt], b);
        return s1.concat(s2);
    } else {
        return [a, b];
    }
}

function findLinearGroups(edge_indices, segments, coordinates) {
    // This routine takes a set of edge_indices and returns a list of 
    // lists, where the edge indices are grouped into lists such that 
    // each list consists of edges that are aligned along a single line
    var toRet = [];
    var edge_remaining = segments.map(x => true);

    for (let ei = 0; ei < edge_indices.length; ei++) {
        let e = edge_indices[ei];
        if (edge_remaining[e]) {
            var current_group = [e];
            var e_pts = segments[e].map(x => coordinates[x]);

            for (let oei = 0; oei < edge_indices.length; oei++) {
                let other_e = edge_indices[oei];

                if ((other_e != e) && (edge_remaining[other_e])) {
                    if (segments[other_e].filter(p => isCollinear(e_pts, coordinates[p])).length > 1) {
                        current_group.push(other_e);
                        edge_remaining[other_e] = false;
                    }
                }
            }
            edge_remaining[e] = false;
            toRet.push(current_group);
        }
    }
    return toRet;
}

function generateInitialRays(R, eis, Li) {
    var all_rays = [];
    var strengths = [];
    for (let i = 0; i < 3; i++) {
        var ei = [0,0,0];
        ei[i] = R;

        var deleted_lattice_points = [[...eis[(i+1)%3]], [...eis[(i+2)%3]]];
        deleted_lattice_points = deleted_lattice_points.concat(Li.filter(x => !x.includes(R)));
        var pts = convexHull(deleted_lattice_points);
	var strPts = pts.map(x => JSON.stringify(x));

        // convexhull includes all lattice points along the sides of the simplex,
        // so we need to remove all of the edge-lattice poitns except the two 
        // closest ones that lie along the side emanating from ei
        var side_point_indices = [0,1,2].map(function(k) {
            return pts.filter(p => p[(i+k)%3] == 0).map(x => strPts.indexOf(JSON.stringify(x)));
	});
        var side1_closest_pt = side_point_indices[1].map(function(pi) {
            var vl = pts[pi].map((p, k) => p - ei[k]);
	    return [veclen(vl), pi];
	}).sort()[0][1];
        var side2_closest_pt = side_point_indices[2].map(function(pi) {
            var vl = pts[pi].map((p, k) => p - ei[k]);
	    return [veclen(vl), pi];
	}).sort()[0][1];

        var side_pts = pts.map(function(p, i) { if (p.includes(0)) {return i;}}).filter(j => j != null);
	var newPts = range(0,pts.length).filter(x => !side_pts.includes(x)).map(y => pts[y]);
	newPts.unshift(pts[side1_closest_pt]);
	newPts.push(pts[side2_closest_pt]);
        var Lis = newPts.map(function(x) {
		return x.map((xp, xpi) => xp - eis[i][xpi]);
	});

        // generate initial rays and strengths for e_i
        var orderedLis = orderedRays(Lis);
        for (let j = 0; j < orderedLis.length; j++) {
            var x = orderedLis[j];
            if (x[1] > 0) {
                all_rays.push([ei, newPts[x[0]]]);
                strengths.push(x[1]);
            }
        }
    }
    return [all_rays, strengths];
}

function identifyRegularTriangles(segments, segment_adjacency_count, coordinates) {
    /*
    * take a set of segments and the number of neighboring triangles for each, 
    * and return a list of 'regular' triangles (triangles whose edges each consist 
    * of r segments, for some r>1
    */
    var ts = [];
    // make a list of the semgents that already have 2 neighbors (i.e. they're done)
    var completed_segments = segment_adjacency_count.map(x => x > 1);
    var interior_segments = segments.map(function(y, i) {if (!(coordinates[y[0]].includes(0) && coordinates[y[1]].includes(0))) {return i;}}).filter(x => x != null);

    // add in all the segments that are on an edge and already have a neighbor
    for (let si = 0; si < segments.length; si++) {
        if ((segment_adjacency_count[si] > 0) && (!interior_segments.includes(si))) {
            completed_segments[si] = true;
        }
    }
    var remaining_segments_index = range(0, segments.length).filter(i => !completed_segments[i]);
    var remaining_segments = remaining_segments_index.map(i => segments[i]);
    var cycles = allCycles(remaining_segments);

    if (cycles.length > 0) {
        for (let ci = 0; ci < cycles.length; ci++) {
            var cycle = cycles[ci];

            var indexed_cycle = cycle.map(x => remaining_segments_index[x]);
            var current_triangle = findLinearGroups(indexed_cycle, segments, coordinates);
            var num_edges_per_side = current_triangle.map(x => x.length);

            // remove cycles that aren't a triangle, cycles that aren't regular triangles, 
            // and the cycle that contains all the outer edges
            if ((indexed_cycle.filter(x => interior_segments.includes(x)).length < 1) || (current_triangle.length != 3) || (Math.max(...num_edges_per_side) - Math.min(...num_edges_per_side) > 0)) {
                //pass;
            } else {
                ts.push(current_triangle.map(ct => ct.map(s => segments[s])));
            }
        }
    }
    return ts;
}

function intersects(seg1, seg2, coordinates) {
    /* 
    * check if two segments intersect at a point interior 
    * to both. i.e. excluding intersection at the endpoints 
    * and intersections of the extended lines passing 
    * through the segments
    */
    var a1 = coordinates[seg1[0]];
    var a2 = coordinates[seg2[0]];
    var b1 = a1.map(function(v, i) {return coordinates[seg1[1]][i] - v;});
    var b2 = a2.map(function(v, i) {return coordinates[seg2[1]][i] - v;});

    if (crossProduct(b1,b2).filter(x => x != 0).length >= b1.length) {
        var t1 = crossProduct([a2[0]-a1[0], a2[1]-a1[1], a2[2]-a1[2]], b2) / crossProduct(b1, b2);
        var t2 = crossProduct([a1[0]-a2[0], a1[1]-a2[1], a1[2]-a2[2]], b1) / crossProduct(b2, b1);

        if ((t1.map(t => t-t1[0]).filter(x => x!=0).length + t2.map(t => t-t2[0]).filter(x => x!=0).length) < 1) {
            if ((t1[0] > 0) && (t1[0] < 1) && (t2[0] > 0) && (t2[0] < 1)) {
                return true;
            }
        }
    }
    return false;
}

function isCollinear(ray, point, tol=0.00001) {
    /*
    * this routine checks if a point lies along the line 
    * defined by the two n-dimensional points: 
    *     ray = [point1, point2]
    */

    // first check if the point is 'equal' to one of the points defining the ray
    if ((point.map(function(p, i) {
	    return (p - ray[0][i])*(p - ray[0][i])}).filter(x => x > tol).length < 1) || (
	 point.map(function(p, i) {
            return (p - ray[1][i])*(p - ray[1][i])}).filter(x => x > tol).length < 1)) {
        return true;
    }
    // otherwise, we'll look at the component-wise slopes
    var numerator = point.map(function(v, i) {return v - ray[0][i];});
    var denominator = ray[1].map(function(v, i) {return v - ray[0][i];});
    var num_nz = numerator.map(v => v*v > 0);
    var den_nz = denominator.map(v => v*v > 0);

    if (num_nz.filter(function(v, i) { return (!(v && den_nz[i]) && (v || den_nz[i]));}).length > 0) {
        return false;
    }
    var frac = numerator.map(function(v, i) {if (num_nz[i]) {return v / denominator[i]}}).filter(y => y != null);

    return (Math.max(...frac) - Math.min(...frac)) < tol;
}

function latticePoints(r=6,a=1,b=2,c=3,nonunit=true) {
    /*
    * Lattice points in the junior simplex lie in the plane 
    * x+y+z = 1, and they are of the form (i/r,j/r,k/r) 
    * where i+j+k=r. 
    * This routine returns the triples with or without scaling by r 
    * (default option for scaling by r is set for nonunit=true)
    */
    var denom = r;
    if (nonunit) {
        denom = 1;
    }
    var points = [[r/denom,0,0],[0,r/denom,0],[0,0,r/denom]];
    for (let i = 0; i < r; i++) {
        let ai = (a*i)%r;
        let bi = (b*i)%r;
        let ci = (c*i)%r;
        if ((ai + bi + ci) == r) {
            points.push(Array.from([ai/denom, bi/denom, ci/denom]));
        }
    }
    return points;
}

function makeTriangulation(edges, boundary_edges=null){
    
    var num_vertices = Math.max(...edges.map(e=>Math.max(e)))+1;
    var nbrs = range(0, num_vertices).map(x => null);
    for (let i = 0; i < edges.length; i++) {
        let e = edges[i];
        if (nbrs[e[0]] != null) {
            if (!nbrs[e[0]].includes(i)) {nbrs[e[0]].push(i);}
        } else {
            nbrs[e[0]] = [i];
        }
        if (nbrs[e[1]] != null) {
            if (!nbrs[e[1]].includes(i)) {nbrs[e[1]].push(i);}
        } else {
            nbrs[e[1]] = [i];
        }
    }

    var triangles = [];
    var edge_to_triangle = edges.map(x=>[]);

    // create set of triples that comprise all possible triangles
    var triples = [];
    var triangle_ctr = 0;
    for (let edge1 = 0; edge1 < edges.length; edge1++) {
        let e = edges[edge1];
        let h = e[0]; let t = e[1];

        // loop over all edges that connect to one endpoint of edge1
        for (let edge2i = 0; edge2i < nbrs[t]; edge2i++) {
            let edge2 = nbrs[t][edge2i];

            if (edge2 != edge1) {
                let p0 = edges[edge2][0]; let p1 = edges[edge2][1];
                let p = (p0 == t) ? p1 : p0;

                // if the neighboring edge connects to one of the neighbors
                // of the second endpoint, then we have a triangle.
                let edge3 = nbrs[p].filter(x => nbrs[h].includes(x));
                if (edge3.length > 0) {
                    edge3 = edge3[0];

                    // check that the edge triple hasn't already been added
                    var st = JSON.stringify([edge1,edge2,edge3].sort());
                    if (!triples.includes(st)) {

                        triples.push(st);

                        triangles.push([edge1,edge2,edge3]);

                        edge_to_triangle[edge1].push(triangle_ctr);
                        edge_to_triangle[edge2].push(triangle_ctr);
                        edge_to_triangle[edge3].push(triangle_ctr);

                        triangle_ctr += 1;
                    }
                }
            }
        }
    }
    /*
    remove the triangles that are extra, and get included because they are subdivisions. 
    That is, triangles like the boundary of the following:
    
                        *
                      / | \
                     /  |  \
                    /  .*.  \
                   / ./   \. \
                  /./       \.\
                 *-------------*
    */
    if (triangles.length > 1) {
        if (boundary_edges != null) {
            if (boundary_edges.length == 3) {
                var boundary_triangle = edge_to_triangle[boundary_edges[0]].filter(function (be0) {
                     return (boundary_edges[1].includes(be0) && boundary_edges[2].includes(be0));
                }).pop();
                if (boundary_triangle != null) {
                    triangles = triangles.splice(boundary_triangle, 1);
                    edge_to_triangle = edge_to_triangle.map(function(x) {
                        return x.map(y => y > extremal_triangle ? y - 1 : y);
                    });
                }
            }
        }
        var ts_to_keep = [];
        for (ti = 0; ti < triangles.length; ti++) {
            let t = triangles[t];
            if (!t.every(e => (edge_to_triangle[e].length != 2))) {
                ts_to_keep.push(ti);
            }
        }
        triangles = ts_to_keep.map(ti => triangles[ti]);
        edge_to_triangle = edge_to_triangle.map(function(x) {
            return x.map(function(y) {
               return ts_to_keep.indexOf(y);
            }).filter(z => z >= 0);
        });
    }
    return [triangles, edge_to_triangle];
}

function matmul(a, b) {
    if (a[0].length != b.length) {
        return null;
    }
    // take transpose of b
    var c = range(0, b[0].length).map(function(i) {
        return b.map(b_row => b_row[i]);
    });
    if (b[0].length < 2) {
	c = [b.map(bj => bj[0])];
    }
    return a.map(function(a_row) {
        return c.map(b_col => dotProduct(a_row, b_col));
    });
}

function nontriangulatedSides(segments, coordinates) {
    /* find out which of the edges contained in the list is not the 
    * side for two triangles. This happens if: 
    * (1) the edge is on the boundary of the domain OR 
    * (2) there are edges missing in the domain, making an incomplete triangulation.
    */
    var all_edge_segments = true;
    var hanging_segments = [];
    var hanging_segment_count = segments.map(x => 0);

    // create quick lookup of all edges emanating from each point
    var node_nbrs = coordinates.map(x => []);
    for (let si = 0; si < segments.length; si++) {
        var s = segments[si];
        node_nbrs[s[0]].push(si);
        node_nbrs[s[1]].push(si);
    }

    // generate all possible triangles (a,b,c) where a,b,c are segments and endpoints agree
    var all_triangles = [];
    for (let s0 = 0; s0 < segments.length; s0++) {
        var s = segments[s0];
        for (let s1i = 0; s1i < node_nbrs[s[0]].length; s1i++) {
            var s1 = node_nbrs[s[0]][s1i];
            if (s1 != s0) {
                for (let s2i = 0; s2i < node_nbrs[s[1]].length; s2i++) {
                    var s2 = node_nbrs[s[1]][s2i];
                    if ((s2 != s1) && (s2 != s0)) {
                        var theSegs = new Set([segments[s1][0], segments[s1][1], segments[s2][0], segments[s2][1]]);
                        if (theSegs.size == 3) {
			    var t = JSON.stringify([s0,s1,s2].sort());
			    if (!all_triangles.includes(t)) {
                                all_triangles.push(t);
			    }
                        }
                    }
                }
            }
        }
    }
    all_triangles = all_triangles.map(function(x) {
        return JSON.parse(x).map(y => parseInt(y));
    });

    var segment_neighbor_count = segments.map(x => 0);
    // now allocate each triangle as a neighbor to its 3 edges
    for (let ti = 0; ti < all_triangles.length; ti++) {
        var t = all_triangles[ti];
        for (let t1i = 0; t1i < t.length; t1i++) {
            var t1 = t[t1i];
            segment_neighbor_count[t1] += 1;
        }
    }

    // tabulate how many segments are missing at least on triangle neighbor
    for (let si = 0; si < segments.length; si++) {
        var s = segments[si];
        if (segment_neighbor_count[si] < 2) {
            hanging_segments.push(s);
            if (s.filter(p => coordinates[p].includes(0)).length < s.length) {
                all_edge_segments = false;
            }
        }
    }
    return [all_edge_segments, hanging_segments, segment_neighbor_count];
}
function nonzeroAvg(a,b,c){
    var cmax = Math.max(...c);
    var idx = c.indexOf(cmax);
    return (a[idx] + b[idx]) / cmax;
}

function orderedRays(L) {
    /*
    * this routine assumes that there exists a Jung-Hirzebruch continued 
    * fraction relation among the rays in the list L, and it returns the
    * order and the 'strength' of each ray in the form of a list with 
    * entries of the form (indexi, ai) to denote the index of L[i] in the 
    * ordered list, and its strength ai. 
    * NOTE: it also assumes that the first and last entry in L 
    * correspond each to a ray at one 'end'
    */
    var l0 = L[0];
    var lk = L[L.length - 1];
    var iLs = L.slice(1, L.length - 1);

    if (iLs.length < 1) {
        return [[0,0],[1,0]];
    } else {
        var angles = argsort(iLs.map(x => veclen(crossProduct(x, l0))));
        var oLs = [l0, ...angles.map(x => iLs[x]), lk]
        var toRet = oLs.slice(0, oLs.length - 2).map(function(oli, i) {
            return [angles[i]+1, nonzeroAvg(oLs[i], oLs[i+2], oLs[i+1])];
        });
        toRet.unshift([0,0]);
        return toRet.concat([L.length - 1, 0]);
    }
}

function orderSegments(segs) {
    var b = segs[0][0];
    var e = segs[0][1];
    var l = [segs[0]];
    var met = segs.map(x => false);
    met[0] = true;

    while(met.filter(x => !x).length > 0) {
        for (let si = 0; si < segs.length; si++) {
            var s = segs[si];
            if (!met[si]) {
                if (s.includes(e)) {
                    var other = s[1 - s.indexOf(e)];
                    l.push([e, other]);
                    e = other;
                    met[si] = true;
                } else {
                    if (s.includes(b)) {
                        var other = s[1 - s.indexOf(b)];
                        l.unshift([other, b]);
                        b = other;
                        met[si] = true;
                    }
                }
            }
        }
    }
    return l;
}

function QPFromTriangulation(t) {
    var ns = [];
    var es = [];
    var fn = [];
    var pt = [];

    if (t != null) {
        var R = Math.max(...t[1].map(x => Math.max(...x)));
        var edges = t[0];
        var coordinates = t[1];
        var extremal_edges = edges.map(function(e, ei) {if (coordinates[e[0]].includes(R) && coordinates[e[1]].includes(R)) {return ei;}}).filter(y => y != null);
        coordinates = rotateSimplexToPlane(coordinates).map(c => [c[0], c[1]]);
        var tri = makeTriangulation(edges, extremal_edges);

        var triangles = tri[0];
        var edge_to_triangle = tri[1];

        var QP_edges = [];
        for (let ti = 0; ti < triangles.length; ti++) {
            var t = triangles[ti];
            for (let i = 1; i < 4; i++) {
                for (let j = -1; j < 2; j += 2) {
                    if (edge_to_triangle[t[i%3]].length > 1 && edge_to_triangle[t[(i+j)%3]].length > 1) {
                        QP_edges.push([t[i%3], t[(i+j)%3]]);
                    }
                }
            }
        }
        for (let i = 0; i < edges.length; i++) {
            var e = edges[i];
            if (edge_to_triangle[i].length > 1) {
                if (curve_type(edges, triangles, edge_to_triangle, coordinates) == [-2,0]) {
                    QP_edges.push([i,i]);
                }
                if (curve_type(edges, triangles, edge_to_triangle, coordinates) == [-3,0]) {
                    QP_edges.push([i,i]);
                    QP_edges.push([i,i]);
                }
            }
        }
        var e_reorder = range(0, edges.length).filter(i => edge_to_triangle[i] > 1);
        QP_edges = QP_edges.map(x => [e_reorder.indexOf(x[0]), e_reorder.indexOf(x[1])]);
        var p = new Set();

        var positions = edges.map(function(e, ei) {
            if (edge_to_triangle[ei].length > 1) {
                return [0.5*(coordinates[e[0]][0] + coordinates[e[1]][0]), 0.5*(coordinates[e[0]][1] + coordinates[e[1]][1])];
            }
        }).filter(y => y != null);

        var xscaling = Math.max(...positions.map(x => x[0])) - Math.min(...positions.map(x => x[0]));
        var yscaling = Math.max(...positions.map(x => x[1])) - Math.min(...positions.map(x => x[1]));

        ns = positions.map(function(p, i) {
    	return {
                "id": i.toString(), "label": i.toString(), 
                "x":(1000.0/xscaling)*p[0], "y":(1000.0/yscaling)*p[1],
    	    };
        });
        es = edges.map(function(e, i) {
            return {
    		"id": i.toString(), "title": "edge "+i.toString(), 
    		"from": e[0].toString(), "to": e[1].toString(), 
    		"arrows": "to"
    	    };
        });
    }
    return JSON.stringify({"nodes": ns, "edges": es, "frozenNodes": fn, "potential": pt});
}

function range(start, stop, step=1) {
    return Array.from({ length: (stop - start) / step}, (_, i) => start + (i * step));
}


function ReidsRecipe(segments, strengths, potential_segments, longest_extension, coordinates) {
    var new_segments = segments.map((s, si) => [-1, s[0], s[1], strengths[si]]);
    var children = segments.map(x => []);

    // make sure each initial segment that meets with an initial segment 
    // of greater strength is marked as "finished" so that it is not 
    // extended further into the domain. 
    var not_finished = strengths.map(x => true);
    for (let i = 0; i < segments.length; i++) {
        var si = segments[i];
        var segs_at = segments.map((sj, j) => sj[1] == si[1] ? j : -1).filter(x => x >=0);
        var strengths_at = segs_at.map(x => strengths[x]);

        for (let ji = 0; ji < segs_at.length; ji++) {
            var j = segs_at[ji];
            not_finished[j] = false;
            var maxAt = Math.max(...strengths_at);
            if ((strengths[j] == maxAt) && (count(strengths_at, maxAt) < 2)) {
                not_finished[j] = true;
            }
        }
    }
    // add in extensions to each remaining non-finished segment
    for (let si = 0; si < segments.length; si++) {
        var s = segments[si];
        if (not_finished[si]) {
            var extensions_for_s = potential_segments.filter(x => x[0] == si).map(
		    function(ps) {return [ps[1], ps[2], 0];
		    });
            if (extensions_for_s.length > 0) {
		var indices_from = [];
		if (extensions_for_s.length > 1) {
                    indices_from = range(0, extensions_for_s.length).map(j => new_segments.length + j);
		}
                indices_from.unshift(si);
                extensions_for_s = extensions_for_s.map(function(e, i) {
                    return [indices_from[i], e[0], e[1], e[2]];
		});
                var added_indices = range(0, extensions_for_s.length).map(j => new_segments.length + j);
                children[si] = [...added_indices];
                for (let e = 0; e < extensions_for_s.length - 1; e++) {
                    children.push(added_indices.slice(e+1));
                }
                children.push([]);
                new_segments.push(...extensions_for_s);
            }
        }
    }
    segments = new_segments;
    var segs_at_pt = coordinates.map(c => []);
    for (let si = 0; si < segments.length; si++) {
	let s = segments[si];
        segs_at_pt[s[2]].push(si);
    }

    var still_updating = true;
    var ctr = 0;
    while((still_updating) && (ctr < longest_extension*segments.length)) {
        ctr += 1;
        still_updating = false;

        // order segments by their strengths
        var current_segs = segments.map((s, si) => [s[s.length - 1], si]).sort().reverse();
        for (let cs = 0; cs < current_segs.length; cs++) {
            var seg_i = current_segs[cs][1];
            var segment = segments[seg_i];
            var seg_i_strength = segments[seg_i][segments[seg_i].length - 1];
            // check if this segment has nonzero strength, and that it has possible extensions
            if ((seg_i_strength >= 0) && (children[seg_i].length > 0)) {
                var segments_at_endpoint = segs_at_pt[segment[2]].filter(s => segments[s][segments[s].length - 1] >= 0);
                var strengths_at_endpoint = segments_at_endpoint.map(s => segments[s][segments[s].length - 1]);
                var first_child = children[seg_i][0];
                if (segment[segment.length - 1] < Math.max(strengths_at_endpoint)) {
                    if (segments[first_child][segments[first_child].length - 1]) {
                        still_updating = true;
                    }
                    for (let ci = 0; ci < children[seg_i].length; ci++) {
                        var c = children[seg_i][ci];
                        segments[c] = [segments[c][0], segments[c][1], segments[c][2], -1];
                    }
                } else {
                    if (count(strengths_at_endpoint, Math.max(strengths_at_endpoint)) < 2) {
                        // add check for intersections here!!!!
                        var child_strength = segment[segment.length - 1] + 1 - segments_at_endpoint.length;
                        if (child_strength != segments[first_child][segments[first_child].length - 1]) {
                            segments[first_child] = [segments[first_child][0], segments[first_child][1], segments[first_child][2], child_strength];
                            still_updating = true;
                        }
                    }
                }
            }
            if (still_updating) {
                break;
            }
        }
    }
    return segments.filter(si => si[si.length - 1] > 0).map(s => [s[1],s[2]]);
}

function rotateSimplexToPlane(coordinates) {
    var n = [1/Math.pow(3,0.5),1/Math.pow(3,0.5),1/Math.pow(3,0.5)];
    var m1 = [
                 [ n[0] / (n[0]*n[0] + n[1]*n[1]), n[1] / (n[0]*n[0] + n[1]*n[1]), 0],
                 [-n[1] / (n[0]*n[0] + n[1]*n[1]), n[0] / (n[0]*n[0] + n[1]*n[1]), 0],
                 [0,0,1]
             ];
    var n1 = matmul(m1, [[n[0]],[n[1]],[n[2]]]);
    var m2 = [
        [n1[2][0], 0, -n1[0][0]],
        [0,1,0],
        [n1[0][0], 0, n1[2][0]]
    ];


    var toRet = coordinates.map(function(c) {
        var c1 = matmul(matmul(m1, m2),[[c[0]], [c[1]], [c[2]]]);
        return [c1[0][0], c1[1][0], c1[2][0]];
    });
    return toRet;
}

function tesselate(triangle, coordinates) {
    var r = triangle[0].length - 1;
    // order the segments along each side so that side1 and side2 emanate from 
    // the same vertex, and so that side3 begins at the end of side1. 
    var side1 = orderSegments(triangle[0]);
    var side2 = orderSegments(triangle[1]);
    var side3 = orderSegments(triangle[2]);
    if ((side3[0][0] == side1[0][0]) || (side3[side3.length - 1][1] == side1[0][0])) {
        var s2 = [...side2];
        side2 = [...side3];
        side3 = s2;
    }
    if (side2[side2.length - 1][1] == side1[0][0]) {
        side2 = side2.reverse().map(y => [y[1],y[0]]);
    }
    if (side3[side3.length - 1][1] == side1[side1.length - 1][1]) {
        side3 = side3.reverse().map(y => [y[1],y[0]]);
    }

    // generate new points
    var new_points = [];
    for (let row = 0; row < r - 1; row++) {
        var num_pts_in_row = r - (row + 1);
        var row_endpoints = [coordinates[side2[row][1]], coordinates[side3[row][1]]];

        if (num_pts_in_row > 1) {
            var dx = (row_endpoints[1][0] - row_endpoints[0][0]) / (num_pts_in_row + 1);
            var dy = (row_endpoints[1][1] - row_endpoints[0][1]) / (num_pts_in_row + 1);
            var dz = (row_endpoints[1][2] - row_endpoints[0][2]) / (num_pts_in_row + 1);
            var interior_points = [...Array(num_pts_in_row).keys()].map(function(i) {
                return [row_endpoints[0][0] + (i+1)*dx, 
                        row_endpoints[0][1] + (i+1)*dy, 
                        row_endpoints[0][2] + (i+1)*dz];
            });
            new_points.push.apply(new_points, interior_points);
        } else {
            if (num_pts_in_row > 0) {
                new_points.push([0.5*(row_endpoints[0][0] + row_endpoints[1][0]), 
                                 0.5*(row_endpoints[0][1] + row_endpoints[1][1]), 
                                 0.5*(row_endpoints[0][2] + row_endpoints[1][2])]);
            }
        }
    }

    // create lookup table o findex points for new_segments and new_points
    let last_entry = side1[side1.length - 1][1];
    var vertex_map = side1.map(s => s[0]).push(last_entry);
    var b = 0;
    for (let row = 0; row < r; row++) {
        vertex_map.push(side2[row][1]);
        if (row < (r - 1)) {
            vertex_map.push.apply(vertex_map, range(b, b + r - (1 + row)).map(x => -(x + 1)));
            b += r - (1 + row);
        }
        vertex_map.push(side3[row][1]);
    }
    vertex_map.push(side2[side2.length - 1][1]);

    // now create segments with appropriately indexed new_pts values
    var new_segments = [];
    var offset = 0;
    for (let row = 0; row < r; row++) {
        new_segments.push(...range(1, r+1-row).map(c => [vertex_map[offset+c], vertex_map[offset+(r+2-row)+c]]));
        new_segments.push(...range(1, r+1-row).map(c => [vertex_map[offset+c], vertex_map[offset+(r+2-row)+c-1]]));
        new_segments.push(...range(1, r+1-row).map(c => [vertex_map[offset+(r+2-row)+c-1], vertex_map[offset+(r+2-row)+c]]));
        offset += r + 2 - row;
    }
    return [new_points, unique(new_segments.map(x => [Math.min(x), Math.max(x)]))];
}

function triangulation(R,a,b,c) {
    // vertices defining the junior simplex (scaled by length R so we deal with integers)
    var eis = [[R,0,0],[0,R,0],[0,0,R]];
    // interior lattice points for junior simplex (also scaled by R)
    var coordinates = latticePoints(R,a,b,c);
    var strCoords = coordinates.map(x => JSON.stringify(x))

    // STEP 1: generate initial rays and their strengths
    var step1output = generateInitialRays(R, eis, coordinates);
    var segments_with_coords = step1output[0];
    var strengths = step1output[1];
    var segments = segments_with_coords.map(r => [strCoords.indexOf(JSON.stringify(r[0])), 
                                                  strCoords.indexOf(JSON.stringify(r[1]))]);

    // STEP 2: claculate possible extensions for each ray
    var potential_segments = [];
    var longest_extension = 0;
    for (let ir = 0; ir < segments_with_coords.length; ir++) {
        var r = segments_with_coords[ir];
        if (strengths[ir] > 0) {
            var points_along_r = coordinates.map(
                function(c, ci) {
                    if (isCollinear(r, c) && !segments[ir].includes(ci)) {
                        return [veclen([c[0]-r[1][0], c[1]-r[1][1], c[2]-r[1][2]]), ci];
	            }
                }).filter(y => y != null);
            points_along_r.unshift([0, segments[ir][1]]);
            potential_segments.push(...points_along_r.map(
                function(v, i) {
                    if (i < points_along_r.length - 1) {
                        return [ir, v[1], points_along_r[i+1][1]];
                    }
                }).filter(y => y != null));
            longest_extension = Math.max(longest_extension, points_along_r.length)
        }
    }

    // STEP 3: find intersections of rays/extensions
    segments = ReidsRecipe(segments, strengths, potential_segments, longest_extension, coordinates);
    // add segments that lie along the sides of the simplex
    for (let i = 0; i < 3; i++) {
        var points_along_side = coordinates.filter(x => x[i] == 0).sort();
        var segments_along_side = points_along_side.map(function(p, pi) {
            if (pi < points_along_side.length - 1) {
                return [strCoords.indexOf(JSON.stringify(p)), strCoords.indexOf(JSON.stringify(points_along_side[pi+1]))];
            }
        }).filter(y => y != null);
        segments.push(...segments_along_side);
        strengths.push(...segments_along_side.map(x => 0));
    }

    // STEP 4: tesselate remaining r-sided triangles into r^2
    var step4outputs = nontriangulatedSides(segments, coordinates);
    var all_edges = step4outputs[0];
    var hanging_segments = step4outputs[1];
    var segment_neighbor_count = step4outputs[2];
    if (!all_edges) {
        var regular_triangles = identifyRegularTriangles(segments, segment_neighbor_count, coordinates);
        for (let ti = 0; ti < regular_triangles.length; ti++) {
            var t = regular_triangles[ti];
            var tesselate_output = tesselate(t, coordinates);
            var extra_points = tesselate_output[0];
            var extra_segments = tesselate_output[1];

            for (let esi = 0; esi < extra_segments.length; sei++) {
                var es = extra_segments[esi];
                var reindexed_segment = [es[0] >= 0 ? es[0] : coordinates.length - (es[0]+1), es[1] >= 0 ? es[1] : coordinates.length - (es[1]+1)];
                segments.push(reindexed_segment);
            }
            if (extra_points.length > 0) {
                coordinates.push(...extra_points);
            }
        }
    }

    return [segments, coordinates];
}

function unique(L) {
    var toRet = new Set(L);
    return Array.from(toRet);
}

function veclen(v) {
    return Math.pow(dotProduct(v,v), 0.5);
}

function viewTriangulation() {
    var localNodes = new vis.DataSet();
    var localEdges = new vis.DataSet();
    if (globalTriangulation != null) {
        var ln = []
        var le = [];

        let rotatedNodes = rotateSimplexToPlane(globalTriangulation[1]);
        var xscaling = Math.max(...rotatedNodes.map(x => x[0])) - Math.min(...rotatedNodes.map(x => x[0]));
        var yscaling = Math.max(...rotatedNodes.map(x => x[1])) - Math.min(...rotatedNodes.map(x => x[1]));
        for (let n = 0; n < rotatedNodes.length; n++) {
    	let xy = rotatedNodes[n];
            ln.push({
    		id: n.toString(),
    		label: n.toString(),
    		x: (1000/xscaling)*parseFloat(xy[0]), 
    		y: (1000/yscaling)*parseFloat(xy[1])
    	});
        }
        for (let e = 0; e < globalTriangulation[0].length; e++) {
    	let s = globalTriangulation[0][e];
            le.push({id: e.toString(), from: s[0].toString(), to: s[1].toString()});
        }

        localNodes.add(ln);
        localEdges.add(le);
    }

    // create a network
    var localContainer = document.getElementById("triangulationView");
    var data = {
        nodes: localNodes,
        edges: localEdges
    };
    var options = {
        nodes: {
            borderWidth:1,
            size:15,
            color: {
                border: '#222222',
                background: 'grey'
            },
 	    physics: {enabled:false},
        },
    };
    localNetwork = new vis.Network(localContainer, data, options);
}
