// triangulation information
var globalBoundaryEdges = [0,1,2,3,4,5];
var globalCoords = [[6,0,0],[0,6,0],[0,0,6], [1,2,3],[2,4,0],[4,2,0],[3,0,3]];
var globalEdges = [[0,5],[5,4],[4,1],[1,2],[2,6],[6,0],[6,5],[5,3],[4,3],[1,3],[2,3],[6,3]];
var globalTriangulation = [globalEdges, globalCoords];

// allows checking if things are the same in euclidean space within this given tolerance
const tolerance = 0.00001;

// information for canvas drawing
var network_nodes, network_edges, triangulate_network;
var network_output_fields = ["id", "from", "to", "edge1","edge2","edge3"] // which data objects to print to screen


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

function allUniqueTriangulations(t, boundary_edges) {
    // return a list of all of the unique triangulations to be obtained by 
    // flipping sequences of edges for the triangulation t. 
    // (the output is represented as JSON stringified lists of edges, 
    // with one such stringifiedlist for each distinct triangulation)

    // order the vertices for each edge so we don't double count edges
    var edges = t[0].map(x => x.sort());
    var coords = t[1];

    // generate initial triangles and edge-to-triangles lists
    var tri = makeTriangulation(edges, coords, boundary_edges);
    var triangles = tri[0].map(x => [...x].sort());
    var edge_to_triangle = tri[1];

    // timeout catch
    var maxRuntime = 10000;
    var beginTime = Date.now();

    // function to convert edges to string of sorted list of sorted lists 
    // (this gives well-posedness when comparing triangulations)
    function sortedString(listOfLists) {
        let tl = listOfLists.map(x => JSON.stringify(x));
        tl.sort();
        return JSON.stringify(tl);
    }

    // recursive function for flipping arbitrarily-long sequences of edges 
    // in the triangulation to obtain "new" triangulations
    function tryFlippingAllEdges(es, ts, e2ts, alreadyMet) {
        if (Date.now() - beginTime > maxRuntime) {
            return;
        }
        for (let e = 0; e < es.length; e++) {
            var flipped = flip(e, es, coords, ts, e2ts, boundary_edges);

            if (flipped[0]) {//check if that edge could be flipped
                var sflipped = sortedString(flipped[1]);
    
                if (alreadyMet.includes(sflipped)) {//check if the current state has already been met
    	        } else {// if it's a new triangulation, add it, and try flipping all edges for this one
    	            alreadyMet.push(sflipped);
                    tryFlippingAllEdges(flipped[1], flipped[2], flipped[3], alreadyMet);
    	        }
    	    }
        }
    }

    var uniqueList = [sortedString(edges)];
    tryFlippingAllEdges(edges, triangles, edge_to_triangle, uniqueList);
    return {triangulations: uniqueList, coordinates: coords, timeout: Date.now() - beginTime > maxRuntime};
}

function argsort(arr) {
    // dsu algorithm for argsort
    let decor = (v, i) => [v, i];
    let undecor = a => a[1];
    let arrcpy = [...arr];
    return arrcpy.map(decor).sort(function(a, b){return a[0]-b[0];}).map(undecor);
}

function getBoundaryEdges (r, edges, coords) {
    console.log(r, edges, coords)
    var gb1 = range(0, edges.length).filter(e => (isCollinear([[r, 0, 0], [0, r, 0]], coords[edges[e][0]]) && isCollinear([[r, 0, 0], [0, r, 0]], coords[edges[e][1]])))
    var gb2 = range(0, edges.length).filter(e => (isCollinear([[0, r, 0], [0, 0, r]], coords[edges[e][0]]) && isCollinear([[0, r, 0], [0, 0, r]], coords[edges[e][1]])))
    var gb3 = range(0, edges.length).filter(e => (isCollinear([[0, 0, r], [r, 0, 0]], coords[edges[e][0]]) && isCollinear([[0, 0, r], [r, 0, 0]], coords[edges[e][1]])))
    return gb1.concat(gb2).concat(gb3)
  }

function callTriangulation() {
    // function that draws initial configuration for triangulate tab of window
    var a = parseFloat(document.getElementById("A_value").value);
    var b = parseFloat(document.getElementById("B_value").value);
    var c = parseFloat(document.getElementById("C_value").value);
    var r = a + b + c;
    var t = triangulation(r,a,b,c);
    globalEdges = t[0];
    globalCoords = t[1];
    globalTriangulation = [t[0], t[1]];

    globalBoundaryEdges = getBoundaryEdges(r, globalEdges, globalCoords)

    viewTriangulation();
}

function canFlip(edge_index, nbr1, nbr2, edges, edge_to_triangle, triangles, coordinates) {
    return curveType(edges, triangles, edge_to_triangle, coordinates, edge_index) < 1;
}

function containsZero(l, tol=tolerance) {
    var squared = l.map(x => x*x);
    return Math.min(...squared) < tol;
}

function convexHull(planar_points) {
    // performs a quickhull search to find which points in the set of 
    // planar_points form the convex hull.
    // This routine assumes that the set of planar_points are points 
    // in 3D, but which lie in some plane. Also assumes that the first 
    // two points in the set of planar_points forms an edge of the 
    // convex hull. The more general version defines the initial starting 
    // points p0 and p1 as any two points with maximal distance between them
    var strVersion = planar_points.map(x => JSON.stringify(x));
    var p0 = planar_points[0];
    var p1 = planar_points[1];
    var s1 = findHull(planar_points.slice(2), p0, p1);
    var s2 = findHull(planar_points.slice(2), p1, p0);
    var points = unique(s1.concat(s2)).map(p => strVersion.indexOf(JSON.stringify(p)));
    return points.map(p => planar_points[p]);
}

function count(l, v) { // return the number of instances of value v in list l
    return l.filter(x => x == v).length;
}

function crossProduct(v1, v2) { // 3d cross product
    return [v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]];
}

function curveType(edges, triangles, edge_to_triangle, coordinates, edge_index) {
    const e = edges[edge_index];
    const e2t = edge_to_triangle[edge_index];
    if (e2t.length < 2) { // check if edge is on the boundary
        return 0;
    }
    var vertex_indices = [];
    for (let ti = 0; ti < 2; ti++) {
	let t = triangles[e2t[ti]];
	for (ei = 0; ei < 3; ei++) {
	    let edgei = edges[t[ei]];
	    if (!vertex_indices.includes(edgei[0])) { vertex_indices.push(edgei[0]);}
	    if (!vertex_indices.includes(edgei[1])) { vertex_indices.push(edgei[1]);}
        }
    }
    let o_verts = vertex_indices.filter(p => !e.includes(p));

    var v1 = coordinates[o_verts[0]];
    var v2 = coordinates[e[0]];
    var v3 = coordinates[o_verts[1]];
    var v4 = coordinates[e[1]];
    
    // solve -Nv2 + Nv4 = v1 + v3 - (v2 + v4)
    var idx = 0;
    for (let i = 0; i < 3; i++) {
        if ( (v4[i] - v2[i])*(v4[i] - v2[i]) > 0 ) {
            idx = i;
            break;
	}
    }
    var rhs = v1[idx] + v3[idx] - (v2[idx] + v4[idx]);
    var N = Math.abs(rhs / (v4[idx] - v2[idx]));
    if (N > 0) {
        N = Math.max(Math.floor(N), Math.floor(1.1/N));
    }
    return N;
}

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

function dotProduct(v1, v2) { // n dimensional dot product (assumes length of v1 = length of v2)
    return v1.map((x, i) => v1[i] * v2[i]).reduce((m, n) => m + n);
}

function findHull(point_set, a, b) { // this is a helper function for the convexhull algorithm
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

function flip(edge_index, edges, coordinates, triangles=[], edge_to_triangle=[], boundary_edges=null, inplace=false) {
    // this routine flips an edge in a given triangulation

    var old_triangles_length = triangles.length;
    if (boundary_edges != null) {
        if (boundary_edges.includes(edge_index)) {
            return [false, edges, triangles, edge_to_triangle];
        }
    }
    var e = edges[edge_index];
    if ((triangles.length < 1) || (edge_to_triangle.length < 1)) {
        var t = makeTriangulation(edges, coordinates, boundary_edges);
        triangles = t[0];
        edge_to_triangle = t[1];
    }
    if (edge_to_triangle[edge_index].length > 1) {
        // find the two triangles that share edge edge_index
        var nbrs = edge_to_triangle[edge_index];
        var nbr1 = nbrs[0]; var nbr2 = nbrs[1];

        if (canFlip(edge_index, nbr1, nbr2, edges, edge_to_triangle, triangles, coordinates)) {
            // get the points p1, p2 'opposite' to edge e
            //    p1-------e[1]
            //    |       /  |
            //    | nbr1 /   |
            //    |     /e   |
            //    |    /     |
            //    |   / nbr2 |
            //    |  /       |
            //   e[0]-------p2
            var nb1_edges = triangles[nbr1].filter(y => y != edge_index);
            var nb2_edges = triangles[nbr2].filter(y => y != edge_index);

            var p1 = edges[nb1_edges[0]].filter(x => !e.includes(x));
            var p2 = edges[nb2_edges[0]].filter(x => !e.includes(x));

            // find the edge indices for e[0]--p2 and e[1]--p1 and other 2 as well.
            var e0p1 = nb1_edges[0];
            var e1p1 = nb1_edges[1];
            var e0p2 = nb2_edges.filter(y => (edges[e0p1].includes(edges[y][0]) || edges[e0p1].includes(edges[y][1])))[0];
            var e1p2 = nb2_edges.filter(y => y != e0p2)[0];

            if (inplace) {
		edges[edge_index] = [Math.min(p1,p2), Math.max(p1,p2)];
		triangles[nbr1] = [e0p2,e0p1,edge_index].sort();
		triangles[nbr2] = [e1p2,e1p1,edge_index].sort();
                edge_to_triangle[e1p1] = edge_to_triangle[e1p1].map(function(y){if(y != nbr1){return y;} else {return nbr2;}});
                edge_to_triangle[e0p2] = edge_to_triangle[e0p2].map(function(y){if(y != nbr2){return y;} else {return nbr1;}});
                return [true, edges, triangles, edge_to_triangle];
	    } else {
                var new_triangles = triangles.map(function(t, ti) {
                    if (ti == nbr1) {
			var a = [e0p2, edge_index, e0p1];
                        return a.sort()
	            } else if (ti == nbr2) {
                        var a = [e1p2, edge_index, e1p1];
                        return a.sort();
	            } else {
                        return [...t].sort();
	            }
	        });
                var new_edges = range(0, edges.length).map(function(edgei) {
                    if ((edgei > edge_index) || edgei < edge_index) {
                        return [...edges[edgei]];
		    } else {
                        return [Math.min(p1,p2), Math.max(p1,p2)];
		    }
		});
		var new_edge_to_triangle = edge_to_triangle.map(function(e2t, e2ti) {
		    if (e2ti == e1p1) {
                        return e2t.map(function(y){if(y != nbr1){return y;} else {return nbr2;}});
		    } else if (e2ti == e0p2) {
                        return e2t.map(function(y){if(y != nbr2){return y;} else {return nbr1;}});
                    } else {
                        return [...e2t];
		    }
		});

                return [true, new_edges, new_triangles, new_edge_to_triangle];
	    }
        }
    }
    return [false, edges, triangles, edge_to_triangle];
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
	}).sort(function(a, b) {return a[0] - b[0];})[0][1];
        var side2_closest_pt = side_point_indices[2].map(function(pi) {
            var vl = pts[pi].map((p, k) => p - ei[k]);
	    return [veclen(vl), pi];
	}).sort(function(a, b) {return a[0] - b[0];})[0][1];

        var side_pts = pts.map(function(p, i) { if (containsZero(p)) {return i;}}).filter(j => j != null);
	var newPts = range(0,pts.length).filter(x => !side_pts.includes(x)).map(y => pts[y]);
	newPts.unshift(pts[side1_closest_pt]);
	newPts.push(pts[side2_closest_pt]);
        var Lis = newPts.map(function(x) {
            return x.map((xp, xpi) => xp - ei[xpi]);
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
    var interior_segments = range(0, segments.length).filter(function(si) {
             let s = segments[si];
             return (!containsZero(coordinates[s[0]])) || (!containsZero(coordinates[s[1]]));
        });

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

function isCollinear(ray, point, tol=tolerance) {
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

function makeTriangulation(edges, coords, boundary_edges=null){
    // this routine takes a set of edges (given as lists of pairs [p1,p2]
    // where p1 and p2 denote the indices for the vertices that make the endpoints of the given edge
    // and it returns a list of all the triangles that can be generated from that configuration 
    // of edges. NOTE: if the set of boundary edges is not provided, then the triangulation might 
    // contain extra triangles (see the comment below)
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

    // create set of triples that comprise all possible triangles
    var triples = [];
    var triangles = [];
    var edge_to_triangle = edges.map(x=>[]);
    var triangle_ctr = 0;
    for (let edge1 = 0; edge1 < edges.length; edge1++) {
        let e = edges[edge1];
        let h = e[0]; let t = e[1];
        // loop over all edges that connect to one endpoint of edge1
        for (let edge2i = 0; edge2i < nbrs[t].length; edge2i++) {
            let edge2 = nbrs[t][edge2i];
            if (edge1 != edge2) {
                let p0 = edges[edge2][0]; let p1 = edges[edge2][1];
		var p = p0;
		if (p0 == t) {
                    p = p1;
		}

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
        var expected = range(0, edges.length).map(e => 2);
        if (boundary_edges != null) {
            if (boundary_edges.length == 3) {
                var boundary_triangle = edge_to_triangle[boundary_edges[0]].filter(function (be0) {
                     return (edge_to_triangle[boundary_edges[1]].includes(be0) && 
			     edge_to_triangle[boundary_edges[2]].includes(be0));
                }).pop();
                if (boundary_triangle != null) {
                    triangles = triangles.filter(function(t, ti) {return ti != boundary_triangle;});
                    edge_to_triangle = edge_to_triangle.map(function(x) {
                        return x.map(function(y) {
                            if (y > boundary_triangle) {
                                return y - 1;
			    } else {
                                if (y < boundary_triangle) {
				    return y;
			        }
			    }
			}).filter(z => z != null);
                    });
                }
            }
            expected = range(0, edges.length).map(function(e) {if (boundary_edges.includes(e)) {return 1;} else {return 2;}});
        }
        var ts_to_keep = [];
        for (let ti = 0; ti < triangles.length; ti++) {
            let t = triangles[ti];
            if (t.every(e => (edge_to_triangle[e].length == expected[e]))) {
                ts_to_keep.push(ti);
            } else {
                let tri_endpoints = unique(t.map(e => edges[e]).flat());
                let other_edges = edges.filter(function(e, ei){return (!t.includes(ei)) && (tri_endpoints.some(p => e.includes(p)))});
                let other_endpoint_coords = unique(other_edges.map(e => [e[0], e[1]]).flat().filter(p => !tri_endpoints.includes(p))).map(c => coords[c]);
                let other_endpoints = other_edges.map(e => [e[0], e[1]]).flat().filter(p => !tri_endpoints.includes(p));

                if (!other_endpoint_coords.some(c => is_interior_to_tri(coords[tri_endpoints[0]], coords[tri_endpoints[1]], coords[tri_endpoints[2]], c))) {
                    ts_to_keep.push(ti);
                }
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

function triangle_area(t1, t2, t3) {
    var x1 = t1[0]; var y1 = t1[1];
    var x2 = t2[0]; var y2 = t2[1];
    var x3 = t3[0]; var y3 = t3[1];
    return Math.abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0);
 }
 
 function is_interior_to_tri(t1, t2, t3, point) {
     //https://www.geeksforgeeks.org/check-whether-a-given-point-lies-inside-a-triangle-or-not/
     /* Calculate area of triangle ABC */
     if ((point[0] > Math.max(t1[0], t2[0], t3[0])) || (point[0] < Math.min(t1[0], t2[0], t3[0]))
      || (point[1] > Math.max(t1[1], t2[1], t3[1])) || (point[1] < Math.min(t1[1], t2[1], t3[1]))) {
         return false;
     }
     var A = triangle_area(t1, t2, t3);
     var A1 = triangle_area(point, t2, t3);
     var A2 = triangle_area(t1, point, t3);
     var A3 = triangle_area(t1, t2, point);
     var B = A1 + A2 + A3;
     return (A == B);
 } 

function matmul(a, b) {
    // matrix multiplication (oh, how I miss numpy...)
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
            if (s.filter(p => containsZero(coordinates[p])).length < s.length) {
                all_edge_segments = false;
            }
        }
    }
    return [all_edge_segments, hanging_segments, segment_neighbor_count];
}

function nonzeroAvg(a,b,c){
    // helper function for orderedRays routine
    var idx = c.indexOf(c.find(x => (x > 0 || x < 0)));
    return (a[idx] + b[idx]) / c[idx];
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
        var angles = iLs.map(x => veclen(crossProduct(x, l0)));
        var indices = argsort(angles);
        var oLs = [l0, ...indices.map(x => iLs[x]), lk];
        var toRet = range(1, oLs.length - 1).map(i => [indices[i-1]+1, nonzeroAvg(oLs[i-1], oLs[i+1], oLs[i])]);
        toRet.unshift([0,0]);
        return toRet.concat([L.length - 1, 0]);
    }
}

function orderSegments(segs) {
    // takes a list of edges (given in pairwise-vertex index format)
    // and returns the list, ordered such that the edges form a conesecutive series
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
    // calculate the quiver from a given triangulation
    var ns = [];
    var es = [];
    var fn = [];
    var pt = [];

    if (t != null) {
        var R = Math.max(...t[1].map(x => Math.max(...x)));
        var edges = t[0];
        var coordinates = t[1];
        var extremal_edges = edges.map(function(e, ei) {
	    e1 = isCollinear([[R,0,0],[0,R,0]], coordinates[e[0]]) && isCollinear([[R,0,0],[0,R,0]], coordinates[e[1]]);
	    e2 = isCollinear([[0,R,0],[0,0,R]], coordinates[e[0]]) && isCollinear([[0,R,0],[0,0,R]], coordinates[e[1]]);
	    e3 = isCollinear([[0,0,R],[R,0,0]], coordinates[e[0]]) && isCollinear([[0,0,R],[R,0,0]], coordinates[e[1]]);
	    if (e1 || e2 || e3) {
                return ei;
	    }
	}).filter(y => y != null);
        var tri = makeTriangulation(edges, coordinates, extremal_edges);

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
            if (edge_to_triangle[i].length > 1) {
		let ct = curveType(edges, triangles, edge_to_triangle, coordinates, i);
		for (let j = 0; j < ct; j++) {
                    QP_edges.push([i, i]);
                }
            }
        }
        var e_reorder = range(0, edges.length).filter(i => edge_to_triangle[i].length > 1);
        QP_edges = QP_edges.map(x => [e_reorder.indexOf(x[0]), e_reorder.indexOf(x[1])]);

        coordinates = rotateSimplexToPlane(coordinates).map(c => [c[0], c[1]]);
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
        es = QP_edges.map(function(e, i) {
            return {
    		"id": i.toString(), "title": "edge "+i.toString(), 
    		"from": e[0].toString(), "to": e[1].toString(), 
    		"arrows": "to"
    	    };
        });
    }
    return JSON.parse(JSON.stringify({"nodes": ns, "edges": es, "frozenNodes": fn, "potential": pt}));
}

function range(start, stop, step=1) { // trying to be as python-ish as I can
    return Array.from({ length: (stop - start) / step}, (_, i) => start + (i * step));
}


function ReidsRecipe(segments, strengths, potential_segments, longest_extension, coordinates) {
    // Reid's recipe helper function
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
                    indices_from = [si, ...range(0, extensions_for_s.length - 1).map(j => new_segments.length + j)];
		}
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
        var current_segs = [...segments].map((s, si) => [s[3], si]);
        current_segs.sort(function(a,b){return a[0] - b[0];}).reverse();

        for (let cs = 0; cs < current_segs.length; cs++) {
            var seg_i = current_segs[cs][1];
            var segment = [...segments[seg_i]];
            var seg_i_strength = segments[seg_i][3];
            // check if this segment has nonzero strength, and that it has possible extensions
            if ((seg_i_strength >= 0) && (children[seg_i].length > 0)) {
                var segments_at_endpoint = segs_at_pt[segment[2]].filter(s => segments[s][3] >= 0);
                var strengths_at_endpoint = segments_at_endpoint.map(s => segments[s][3]);
                var first_child = children[seg_i][0];
                if (segment[3] < Math.max(...strengths_at_endpoint)) {
                    if (segments[first_child][3] >= 0) {
                        still_updating = true;
                    }
                    for (let ci = 0; ci < children[seg_i].length; ci++) {
                        var c = children[seg_i][ci];
                        segments[c] = [segments[c][0], segments[c][1], segments[c][2], -1];
                    }
                } else {
                    if (count(strengths_at_endpoint, Math.max(...strengths_at_endpoint)) < 2) {
                        // add check for intersections here!!!!
                        not_intersects = true;
                        for (let sj = 0; sj < segments.legnth; sj++) {
                            let s = segments[sj];
                            if ((sj != first_child) && (s[3] > 0)) {
                                if (intersects([s[1], s[2]], [segments[first_child][1], segments[first_child][2]], coordinates)) {
                                    not_intersects = false; break;
                                }
                            }
                        }
                        if (not_intersects) {
                            var child_strength = segment[3] + 1 - segments_at_endpoint.length;
                            if (child_strength != segments[first_child][3]) {
                                segments[first_child] = [segments[first_child][0], segments[first_child][1], segments[first_child][2], child_strength];
                                still_updating = true;
                            }
                        } else {
                            if (segments[first_child][-1] >= 0) {
                                still_updating = true;
                            }
                            for (let ci = 0; ci < children[seg_i].length; ci++) {
                                let c = children[seg_i][ci];
                                segments[c] = [segments[c][0], segments[c][1], segments[c][2], -1];
                            }
                        }
                    }
                }
            }
            if (still_updating) {
                break;
            }
        }
    }
    return segments.filter(function(si) {return si[3] > 0;}).map(s => [s[1], s[2]]);
}

function resolveNetworkFlip(n, p) {
    // this routine flips the clicked edge on the canvas, and updates globals accordingly
    var e = p.edges[0].toString();
    var tri = makeTriangulation(globalTriangulation[0], globalCoords, globalBoundaryEdges);
    var opt = flip(parseInt(e), globalTriangulation[0], globalTriangulation[1], tri[0], tri[1], globalBoundaryEdges);
    if (opt[0]) {
	var es = JSON.parse(JSON.stringify(opt[1])).map(
            function(x) {
	        return x.map(y => parseInt(y));
	    });
	globalEdges = es;
	globalTriangulation = [globalEdges, globalCoords];

    updateNetworkTriFromGlobal();

    document.getElementById("tri-edges").innerText = JSON.stringify(network_edges.get(),output_fields, 4);
    var tri = opt[2].map(function(t) {
            return {
                edge1: "index: "+t[0].toString() + ", "+JSON.stringify(globalEdges[t[0]]), 
                edge2: "index: "+t[1].toString() + ", "+JSON.stringify(globalEdges[t[1]]), 
                edge3: "index: "+t[2].toString() + ", "+JSON.stringify(globalEdges[t[2]])
            };
    });
    document.getElementById("tri-triangles").innerText = JSON.stringify(tri,output_fields, 4);
    updateGlobalTriFromNetwork();
    }
}

function rotateSimplexToPlane(coordinates) {
    // helper function for viewing the canvas
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

function showTriExchangeNumber() {
    const output = document.getElementById('tri-exchange-number-output')
    try {
        const result = allUniqueTriangulations(globalTriangulation, globalBoundaryEdges);
        if (result.timeout) {
            output.textContent = "Timed out"
        } else {
            output.textContent = result.triangulations.length;
        }
    } catch(e) {
        console.error(e);
        output.textContent = "Error"
    }
}

function tesselate(triangle, coordinates) {
    // this routine takes a regular triangle and subdivides it into unit-segment triangles. 
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
        side2.reverse();
        side2 = side2.map(y => [y[1],y[0]]);
    }
    if (side3[side3.length - 1][1] == side1[side1.length - 1][1]) {
        side3.reverse();
        side3 = side3.map(y => [y[1],y[0]]);
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
            var interior_points = range(0, num_pts_in_row).map(function(i) {
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
    var vertex_map = side1.map(s => s[0]).concat(last_entry);
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
    return [new_points, unique(new_segments.map(x => [Math.min(...x), Math.max(...x)]))];
}

function triangulation(R,a,b,c) {
    /* this routine generates the G-Hilb triangulation for the group 1/R(a,b,c) */
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
            points_along_r.sort(function(a, b) {return a[0] - b[0]});
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
        var points_along_side = coordinates.filter(x => x[i] == 0).sort(function(a, b) {return a[(i+1)%3]-b[(i+1)%3]});
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

            for (let esi = 0; esi < extra_segments.length; esi++) {
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

function updateNetworkFromGlobals() {
    // update the canvas to mirror the global edge/coordinates data
    network_nodes.clear();
    network_edges.clear();

    let rotatedNodes = rotateSimplexToPlane(globalTriangulation[1]);
    let xscaling = Math.max(...rotatedNodes.map(x => x[0])) - Math.min(...rotatedNodes.map(x => x[0]));
    let yscaling = Math.max(...rotatedNodes.map(x => x[1])) - Math.min(...rotatedNodes.map(x => x[1]));
    for (let n = 0; n < rotatedNodes.length; n++) {
        let xy = rotatedNodes[n];
        network_nodes.add({
            id: n.toString(),
            label: n.toString(),
            x: -50 + (300/xscaling)*parseFloat(xy[0]),
            y: -50 + (300/yscaling)*parseFloat(xy[1])
	});
    }
    for (let e = 0; e < globalTriangulation[0].length; e++) {
        let s = globalTriangulation[0][e];
        network_edges.add({id: e.toString(), from: s[0].toString(), to: s[1].toString()});
    }
}

function veclen(v) { // returns the length of a vector v
    return Math.pow(dotProduct(v,v), 0.5);
}

function viewTriangulation() {
    network_nodes = new vis.DataSet();
    network_edges = new vis.DataSet();

    if (globalTriangulation != null) {
        var ln = []
        var le = [];

        let rotatedNodes = rotateSimplexToPlane(globalTriangulation[1]);
	let xscaling = Math.max(...rotatedNodes.map(x => x[0])) - Math.min(...rotatedNodes.map(x => x[0]));
	let yscaling = Math.max(...rotatedNodes.map(x => x[1])) - Math.min(...rotatedNodes.map(x => x[1]));
        for (let n = 0; n < rotatedNodes.length; n++) {
    	    let xy = rotatedNodes[n];
            ln.push({
    		id: n.toString(),
    		label: n.toString(),
    		x: -50 + (300/xscaling)*parseFloat(xy[0]),
    		y: -50 + (300/yscaling)*parseFloat(xy[1])
    	    });
        }
        for (let e = 0; e < globalTriangulation[0].length; e++) {
    	    let s = globalTriangulation[0][e];
            le.push({id: e.toString(), from: s[0].toString(), to: s[1].toString()});
        }

        network_nodes.add(ln);
        network_edges.add(le);
    }

    // create a network
    var localContainer = document.getElementById("triangulationView");
    var data = {
        nodes: network_nodes,
        edges: network_edges
    };
    var options = {
        nodes: {
            borderWidth:1,
            size:10,
            color: {
                border: '#222222',
                background: 'grey'
            },
 	    physics: {enabled:false},
        },
        edges: {
            "smooth": {
                "type": "continuous",
                "forceDirection": "none",
                "roundness": 0
            },
	},
    };
    triangulate_network = new vis.Network(localContainer, data, options);
    triangulate_network.on('click',function(params){
        resolveNetworkFlip(triangulate_network, params);
    });

    document.getElementById("tri-edges").innerText = JSON.stringify(network_edges.get(), network_output_fields, 4);
    var opt = makeTriangulation(globalEdges, globalCoords, globalBoundaryEdges);
    var tri = opt[0].map(function(t) {
	    return {
	        edge1: "index: "+t[0].toString() + ", "+JSON.stringify(globalEdges[t[0]]), 
                edge2: "index: "+t[1].toString() + ", "+JSON.stringify(globalEdges[t[1]]), 
                edge3: "index: "+t[2].toString() + ", "+JSON.stringify(globalEdges[t[2]])
	    };
    });
    document.getElementById("tri-triangles").innerText = JSON.stringify(tri, network_output_fields, 4);
}
