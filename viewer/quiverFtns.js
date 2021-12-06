var nodes, edges, network;
var frozen_nodes;
var potential;
var output_fields = ["id", "from", "to", "coef"] // which data objects to print to screen

function updateInstructions() {
    let click_mode = document.getElementById("edit-quiver-type").value;
    var textOutput = "Click on canvas to add a node";

    if (click_mode == "add-edge") {
	textOutput = "Select a pair of nodes to add an edge between them: control + click to select two nodes";
    } else if (click_mode == "add-loop") {
	textOutput = "Click on a node to add a loop";
    } else if (click_mode == "remove-edge") {
	textOutput = "Click on an edge to remove it";
    } else if (click_mode == "add-node") {
	textOutput = "Click on the canvas to create a node";
    } else {
        textOutput = "Select a node to perform action";
    }
    document.getElementById("instructions").innerText = textOutput;
}

function resolve_click_event(n, p) {
    var click_mode = "NONE";
    try {
        click_mode = document.getElementById("edit-quiver-type").value;
    } catch (err) {
        alert(err)
    }

    if (click_mode == "add-node") {
        if((p.nodes.length == 0) && (p.edges.length == 0)) { // make sure we're not superimposing nodes
	    var nv = getUniqueNodeId();
            nodes.add({id:nv, label:nv, x:p.pointer.canvas.x, y:p.pointer.canvas.y});
        }
    } else if (click_mode == "remove-node") {
        if (p.nodes.length > 0) {
            nodes.remove({id: p.nodes[0].toString()});
        }
    } else if (click_mode == "add-edge") {
        var ne = getUniqueEdgeId();
	if (p.nodes.length > 1) {
	    let n1 = p.nodes[0].toString();
	    let n2 = p.nodes[1].toString();
            edges.add({id: ne, from: n1, to: n2, arrows: "to"});
        }
    } else if (click_mode == "remove-edge") {
        if (p.edges.length > 0) {
            edges.remove({id: p.edges[0].toString()});
        }
    } else if (click_mode == "add-loop") {
        var ne = getUniqueEdgeId();
        if (p.nodes.length > 0) {
	    let v = p.nodes[0].toString();
            edges.add({id: ne, from: v, to: v});
        }
    } else if (click_mode == "freeze-node") {
        if (p.nodes.length > 0) {
            frozen_nodes.add({id: p.nodes[0].toString()});
        }
    } else if (click_mode == "unfreeze-node") {
        if (p.nodes.length > 0) {
            frozen_nodes.remove({id: p.nodes[0].toString()});
        }
    } else if (click_mode == "mutate") {
        if (p.nodes.length > 0) {
            let current_node = p.nodes[0].toString();
            var mQP = makeQP(edges, nodes, frozen_nodes, potential);

            mQP = mutateQP(current_node, mQP);
            var outputs = [];
            for (let i = 0; i < mQP.edges.length; i++) {
                outputs.push({id: i.toString(), to: mQP.edges[i][1].toString(), from: mQP.edges[i][0].toString(), arrows: "to"});
            }
            edges.clear();
            edges.add(outputs);

	    potential.clear();
            potential.add(mQP.potential.filter(x => (x[1] != ",")).map(function(x) {
                return {
                    id: x[1],
                    coef: x[0].toString(),
                };
            }));
        }
    }
}

function clearPotential() {
    potential.clear();
}

function clearQP() {
    edges.clear();
    frozen_nodes.clear();
    nodes.clear();
    potential.clear();
}

function updatePotential() {
    var t = document.getElementById("term-input").value;
    var c1 = parseInt(document.getElementById("coefficient-input").value);
    var c2 = 0;
    try {
        c2 = potential.get(t);
    } catch(err) {
    }
    if (c2 == null) {
        potential.add({id: t, coef: c1.toString()});
    } else {
	let c3 = c1 + parseInt(c2.coef);
        potential.update({id: t, coef: c3.toString()});
    }
}

function generateRandomPotential() {
    var np = randomPotential(nodes, edges);

    potential.clear();
    potential.add(np.filter(x => (x[1] != ",")).map(function(x) {
    return {
        id: x[1],
        coef: x[0].toString(),
        };
    }));
}

function getUniqueEdgeId() { /* create a string value that is not currently in ids for edges */
    var ne = edges.getIds();
    ne = ne.map(function(x) {return parseInt(x);});
    return ne.reduce(function(a, b) {return Math.max(a, b) + 1;}, 0).toString();
}

function getUniqueNodeId() { /* create a string value that is not currently in ids for nodes */
    var nv = nodes.getIds();
    nv = nv.map(function(x) {return parseInt(x);});
    return nv.reduce(function(a, b) {return Math.max(a, b) + 1;}, 0).toString();
}

function draw() {
    // create an array with nodes
    nodes = new vis.DataSet();
    nodes.on("*", function () {
        document.getElementById("nodes").innerText = JSON.stringify(
            nodes.get(),
            output_fields,
            4
        );
    });
    nodes.add([
        { id: "0", label: "0", x:-100, y: 0},
        { id: "1", label: "1", x: 0, y: 100},
        { id: "2", label: "2", x: 100, y: 0},
        { id: "3", label: "3", x: 100, y:-100},
        { id: "4", label: "4", x: 0, y:-200},
        { id: "5", label: "5", x:-100, y:-100},
    ]);

    // create an array with edges
    edges = new vis.DataSet();
    edges.on("*", function () {
        document.getElementById("edges").innerText = JSON.stringify(
            edges.get(),
            output_fields,
            4
        );
    });
    edges.add([
        { id: "0", from: "0", to: "1", arrows: "to"},
        { id: "1", from: "1", to: "0", arrows: "to"},
        { id: "2", from: "1", to: "2", arrows: "to"},
        { id: "3", from: "2", to: "1", arrows: "to"},
        { id: "4", from: "0", to: "2", arrows: "to"},
        { id: "5", from: "2", to: "0", arrows: "to"},
        { id: "6", from: "2", to: "3", arrows: "to"},
        { id: "7", from: "3", to: "2", arrows: "to"},
        { id: "8", from: "3", to: "3", arrows: "to"},
        { id: "9", from: "3", to: "4", arrows: "to"},
        { id: "10", from: "4", to: "3", arrows: "to"},
        { id: "11", from: "4", to: "4", arrows: "to"},
        { id: "12", from: "4", to: "5", arrows: "to"},
        { id: "13", from: "5", to: "4", arrows: "to"},
        { id: "14", from: "5", to: "5", arrows: "to"},
        { id: "15", from: "5", to: "5", arrows: "to"},
        { id: "16", from: "5", to: "0", arrows: "to"},
        { id: "17", from: "0", to: "5", arrows: "to"},
    ]);
    frozen_nodes = new vis.DataSet();
    frozen_nodes.on("*", function () {
    document.getElementById("frozen_nodes").innerText = JSON.stringify(
          frozen_nodes.get(),
          output_fields,
	  4);
    });
    potential = new vis.DataSet();
    potential.on("*", function () {
    document.getElementById("potential").innerText = JSON.stringify(
          potential.get(),
          output_fields,
	  4);
    });
    potential.add([
        { id: "0,2,5", coef: "1"},
        { id: "1,4,3", coef: "1"},
        { id: "8,9,10", coef: "1"},
        { id: "2,6,7,3", coef: "1"},
        { id: "0,1,17,16", coef: "1"},
        { id: "6,8,7", coef: "1"},
    ]);

    // create a network
    var container = document.getElementById("mynetwork");
    var data = {
        nodes: nodes,
        edges: edges,
    };
    var options = {
        nodes: {
            borderWidth:1,
            size:45,
            color: {
                border: '#222222',
                background: 'grey'
            },
            font:{color:'black',
            size: 11,
            face :'arial',
            },
 	   physics: {enabled:false},
        },
        edges: {
            arrows: {
                to:{enabled: true},
            },
            color: {
                color:'#848484',
                highlight:'#848484',
                hover: '#848484',
            },
            font: {
                color: '#343434',
                size: 11, // px
                face: 'arial',
                background: 'none',
                strokeWidth: 5, // px
                strokeColor: '#ffffff',
                align:'vertical'
            },
 	   //physics: {enabled:true},
        },
	interaction: { multiselect: true},
	navigation: true,
     };
     network = new vis.Network(container, data, options);
     //network.setOptions({physics:{enabled:false}});
     network.on('click',function(params){
	 resolve_click_event(network, params);
     });
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*        Begin QP manipulation functions (no global variables)          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
function allThreeCycles (QP) {
    var triples = [];
    for (let v1 in QP.nodes) {
        let arrowsOut1 = QP.arrowsWithTail[v1];

        for (let e1ii = 0; e1ii < arrowsOut1.length; e1ii++) {
            let e1i = arrowsOut1[e1ii];
            let e1 = QP.edges[e1i];
            let v2 = e1[1];
            let arrowsOut2 = QP.arrowsWithTail[v2];

            for (let e2ii = 0; e2ii < arrowsOut2.length; e2ii++) {
                let e2i = arrowsOut2[e2ii];

                if (e2i != e1i) {
                    let e2 = QP.edges[e2i];
                    let v3 = e2[1];
                    let arrowsOut3 = QP.arrowsWithTail[v3];

                    for (let e3ii = 0; e3ii < arrowsOut3.length; e3ii++) {
                        let e3i = arrowsOut3[e3ii];

                        if ((e3i != e1i) && (e3i != e2i)) {
                            let e3 = QP.edges[e3i];
                            let v4 = e3[1];

                            if (v4 == v1) {
                                let triple = cycleOrder([e1i, e2i, e3i]).toString();
                                if (!triples.includes(triple)) {
                                    triples.push([e1i,e2i,e3i]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return triples;
}

function arrayEquals(a, b) {
    return JSON.stringify(a) == JSON.stringify(b);
}

function cycleOrder(cycle) {
    // order a list of integers with minimal element first 
    // (note: need to fix this for multiple instances of min element)
    let thisCycle = cycle.filter(y => y != null).map(x => parseInt(x));
    let minVal = Math.min(...thisCycle);
    let minIdx = thisCycle.indexOf(minVal);
    return thisCycle.slice(minIdx).concat(thisCycle.slice(0, minIdx));
}

function deepCopy(A) {
    return JSON.parse(JSON.stringify(A));
}

function edgesOutOf(vertex, edge_list) {
    return Array.from(Array(edge_list.length).keys()).map(x => parseInt(x)).filter(x => edge_list[x][0] == vertex);
}

function findCycleDFS(begin_vertex, at_vertex, edge_list, edges_so_far, min_length=0) {
    let edgesOut = edgesOutOf(at_vertex, edge_list).filter(y => !edges_so_far.includes(y));

    for (let ei = 0; ei < edgesOut.length; ei++) {
        let e = edge_list[edgesOut[ei]];
	let esMet = deepCopy(edges_so_far).map(y => parseInt(y));
        esMet.push(edgesOut[ei]);
	let currentAt = e[1];
	              
        if ((begin_vertex == currentAt) && (esMet.length >= min_length)) {
            return [true, esMet];
	}
	return findCycleDFS(begin_vertex, currentAt, edge_list, esMet, min_length);
    }
    return [false, esMet];
}

function findDependencies(currentItem, met, lookups) {
    /* find the items that currentItem depends on in the lookup table 
     * (this is recursively done, so the tree of dependencies is listed exhaustively)
     * stopping criterion: no new dependencies found. 
     * circularity is dealt with by only adding new dependants */
    if ( met.includes(currentItem) ) {
        return (true, met);
    } else {
        var newMet = met.concat(currentItem);
        for (let i = 0; i < lookups[currentItem].length; i++) {
	    let item = lookups[currentItem][i];
	    return findDependencies(item, newMet, lookups)
        }
        return findDependencies(false, met);
    }
}

function makeQP(es, ns, fn, p, inputType="fromVisDataSet") {
    // create graph info (arrow to node structures)
    var awh = Array.from(Array(ns.length), x => []);
    var awt = Array.from(Array(ns.length), x => []);
    var la = Array.from(Array(ns.length), x => []);
    var fns = [];
    var theseNodes = [];
    var theseEdges = [];
    var thisPotential = [];

    if (inputType == "fromVisDataSet") {
        fns = fn.getIds().map(x => parseInt(x));
        theseNodes = ns.getIds().map(x => parseInt(x));
        theseEdges = es.getIds().map(x => [parseInt(es.get(x).from), parseInt(es.get(x).to)]);
        // make sure that the edges in each cycle of the potential are now listed by their index,
        // rather than id, since adding/subtracting edges can leave gaps in the id#s 
	var edgeIDMap = es.getIds();
        thisPotential = p.getIds().map(x => [parseInt(p.get(x).coef), x.split(",").map(y => edgeIDMap.indexOf(y)).toString()]);
    } else {
        fns = Array.from(fns, x => parseInt(x));
        theseNodes = Array.from(ns, x => parseInt(x));
        theseEdges = deepCopy(es.filter(x => (x != null))).map(x => [parseInt(x[0]), parseInt(x[1])]);
        thisPotential = p;
    }

    for (let ei = 0; ei < theseEdges.length; ei++) {
        let e = theseEdges[ei];
        if (e[0] == e[1]) {
            la[theseNodes.indexOf(e[0])].push(ei);
        } else {
            awh[theseNodes.indexOf(e[1])].push(ei);
            awt[theseNodes.indexOf(e[0])].push(ei);
        }
    }
    // can_mutate list
    var cm = theseNodes.map(
        function(v) {
            return ((la[v].length < 1) && !(v in fns));
        });
    return {
        arrowsWithHead: awh,
        arrowsWithTail: awt,
        canMutate: cm,
        edges: theseEdges,
        frozenNodes: deepCopy(fns),
        loopsAt: la,
        nodes: theseNodes,
        potential: thisPotential
    }
}

function mutateQP(vertex, QP) {
    const v = parseInt(vertex);
    if (QP.canMutate[v]) {
        // reverse the arrows incident to vertex
        var savedEdges = deepCopy(QP.edges);
        for (let ii in QP.arrowsWithHead[v]) {
            let i = QP.arrowsWithHead[v][ii];
            let e = savedEdges[i];
            savedEdges[parseInt(i)] = [e[1],e[0]];
        }
        for (let ii in QP.arrowsWithTail[v]) {
            let i = QP.arrowsWithTail[v][ii];
            let e = savedEdges[i];
            savedEdges[parseInt(i)] = [e[1],e[0]];
        }

        var delta = [];
        var edgeCtr = QP.edges.length;
        var shortcuts = [];
        // now add 'shortcuts'
        for (let ei1i in QP.arrowsWithHead[v]) {
            let ei1 = QP.arrowsWithHead[v][ei1i];
            var i = QP.edges[parseInt(ei1)][0];
            for (let ei2i in QP.arrowsWithTail[v]) {
                let ei2 = QP.arrowsWithTail[v][ei2i];
                var j = QP.edges[ei2][1];

                shortcuts.push([ei1,ei2]);
                savedEdges.push([i,j]);
                delta.push([1, ei2.toString()+","+ei1.toString()+","+edgeCtr.toString()]);
                edgeCtr++;
            }
        }

        var wPrime = [];
        // update the potential
        for (let mci = 0; mci < QP.potential.length; mci++) {
            var coef = QP.potential[mci][0];
            var monoid = QP.potential[mci][1].split(',');
            let ml = monoid.length;
            var m = "";
            var foundMatch = false;
            for (let i = 0; i < monoid.length; i++) {
                m1 = parseInt(monoid[i]);
                m2 = parseInt(monoid[(i+1)%ml]);
                m0 = parseInt(monoid[(i+ml-1)%ml]);

                const isIn = shortcuts.findIndex(e => arrayEquals(e, [m1, m2]));
                const wasIn = shortcuts.findIndex(e => arrayEquals(e, [m0, m1]));
                if ((i > 0) || (wasIn < 0)) {
                    if((isIn >= 0) && !foundMatch) {
                        var val = QP.edges.length + parseInt(isIn);
                        m = m + val.toString()+",";
                        foundMatch = true;
                    } else {
                        if (!foundMatch) {
                            m = m + m1.toString()+",";
                        }
                        foundMatch = false;
                    }
                }
            }
            wPrime.push([coef, m.slice(0, -1)]);
        }

        wPrime.push(...delta);
        // reduce the resulting quiver
        return reduce(makeQP(savedEdges, QP.nodes, QP.frozenNodes, wPrime, inputType="fromQP"));
    } else {
        return makeQP(QP.edges, QP.nodes, QP.frozenNodes, QP.potential, inputType="fromQP");
    }
}

function pathDerivative(thisPotential, edgeIndex) {
    return thisPotential.map(
	function(x) {
            const inTerm = x[1].indexOf(edgeIndex.toString());
            if (inTerm >= 0) {
                const partial = x[1].split(',').map(
                    function(y) {
			if (y != edgeIndex.toString()) {
			    return parseInt(y);
			}
		    });
                return [x[0], partial];
	    }
	}).filter(y => (y != null));
}

function randInt(range) {
    let sign = Math.random() > 0.5 ? 1 : -1;
    return sign * Math.floor(Math.random() * range);
}

function randomPotential(ns, es, coefficient_range=100) {
    theseNodes = ns.getIds().map(x => parseInt(x));
    theseEdges = es.getIds().map(x => [parseInt(es.get(x).from), parseInt(es.get(x).to)]);
    let thisQP = makeQP(theseEdges, theseNodes, [0], [[0,"0"]], "fromThing");
    let allTriples = allThreeCycles(thisQP);

    var currentlyMetNodes = [];
    var metAllNodes = false;
    var cycles = [];
    var metCycles = [];
    do {
        let randomCycle = randInt(allTriples.length);
        randomCycle = randomCycle >= 0 ? randomCycle : -randomCycle;
        let theCycle = allTriples[randomCycle];

        if (!metCycles.includes(randomCycle)) {
            for (let tci = 0; tci < 3; tci++) {
                let e = theseEdges[theCycle[tci]];
                if (!currentlyMetNodes.includes(e[0])) {
                    currentlyMetNodes.push(e[0]);
                }
                if (!currentlyMetNodes.includes(e[1])) {
                    currentlyMetNodes.push(e[1]);
                }
            }
            metCycles.push(randomCycle);
            cycles.push(theCycle.toString());
        }

        if ((currentlyMetNodes.length >= ns.length) || (metCycles.length >= allTriples.length)){
            metAllNodes = true;
        }
    } while(!metAllNodes);

    return cycles.map(x => [randInt(coefficient_range), x]);
}

/*
function randomPotential(ns, es, coefficient_range=100) {
    theseNodes = ns.getIds().map(x => parseInt(x));
    theseEdges = es.getIds().map(x => [parseInt(es.get(x).from), parseInt(es.get(x).to)]);
    var cycles = [];
    var currentPath = [];

    // find n-cycles containing each node, where n is the
    // smallest number >=3 that admits a cycle
    for (let i = 0; i < theseNodes.length; i++) {
        var output = findCycleDFS(i, i, theseEdges, currentPath, min_lenth=3);
	if (output[0]) {
            let cycle = cycleOrder(output[1]).toString();
            if (!cycles.includes(cycle)) {
                cycles.push(cycle);
            }
	}
    }
    return cycles.map(x => [randInt(coefficient_range), x]);
}
*/

function reduce(QP) {
    // remove extraneous commas from potential
    var thePotential = QP.potential.map(
        function(x) {
            var y = x[1];
            if (y[0] == ",") { y = y.slice(1); }
            if (y[-1] == ",") { y = y.slice(0,-1); }
            return [x[0], y];
        });
    // extract list of unique entries all edges that occur in quadratic terms in the potential
    var squareTerms = thePotential.filter(x => x[1].split(",").length == 2).map(y => y[1]);
    var squareCoefs = thePotential.filter(x => x[1].split(",").length == 2).map(y => y[0]);

    // if there are enough quadratic terms to cancel them
    if (squareTerms.length > 0) {
        var edgesToRemove = new Array();
        var reduceDict = [...Array(QP.edges.length).keys()].map(x => [[1, [parseInt(x)]]]);
        for (let ti = 0; ti < squareTerms.length; ti++) {
            let t = squareTerms[ti];
            let e1 = parseInt(t.slice(0, t.indexOf(',')));
            let e2 = parseInt(t.slice(t.indexOf(',')+1));
            let c = parseInt(squareCoefs[ti]);
            if (!edgesToRemove.includes(e1)) { edgesToRemove.push(e1); }
            if (!edgesToRemove.includes(e2)) { edgesToRemove.push(e2); }

            reduceDict[e1] = pathDerivative(thePotential, parseInt(e2)).map(
                function(x) {
                    if ((x[1].filter(y => y!=null).length > 1) || !x[1].includes(e1)) {
                        return [-parseInt(x[0])/c, x[1]];
                    }
                }).filter(y => (y != null));
            reduceDict[e2] = pathDerivative(thePotential, parseInt(e1)).map(
                function(x) {
                    if ((x[1].filter(y => y!=null).length > 1) || !x[1].includes(e2)) {
                        return [-parseInt(x[0])/c, x[1]];
                    }
                }).filter(y => (y != null));
        }
	// check if there is a circular dependency between any of the quadratic terms 
        // (this makes it so that one edge at least cannot be deleted in the reduction process)
        if (edgesToRemove.length%2 > 0) {
            let lookups = edgesToRemove.map(function(x) {
                var toRet = [];
                for (let yi = 0; yi < reduceDict[x]; yi++) {
                    let y = reduceDict[yi][1];
                    toRet.push(y.filter(z => edgesToRemove.includes(z)));
                }
                return toRet;
	    });
	    var loopStarts = [];
	    var added = [];
            for (let ei = 0; ei < edgesToRemove.length; ei++) {
                let e = edgesToRemove[ei];
	        var opt = [];
	        let outcomes = findDependencies(ei, opt, lookups);
                if (outcomes[0]) {
	            var o1 = outcomes[1];
		    o1.sort();
		    let c1 = o1[0];
		    o1 = o1.toString();
                    if (!added.includes(o1)) {
                        loopStarts.push(parseInt(c1));
			added.push(o1);
	            }
	        }
	    }
	    if (loops.length > 0) {
                edgesToRemove = edgesToRemove.filter(x => !loopStarts.includes(x));

	        for (let i = 0; i < loopStarts.length; i++) {
                    reduceDict[loopStarts[i]] = [[1, [loopStarts[i]]]];
	        }
	    }
	}

        // update lookup table of replacements for each edge to be removed so that
        // each replacement term only contains edges that are not in edgesToRemove
        for (let etri = 0; etri < edgesToRemove.length; etri++) {
            let e = parseInt(edgesToRemove[etri]);
            var termsForE = deepCopy(reduceDict[e]);
            var foundReplacement = true;
            var ctr = 0;

            do {
		// stopping criteria
                foundReplacement = true; ctr++;

                //placeholder for holding non-edgesToRemove lookup values for edge e
                var altTermsForE = []
                for (let cti = 0; cti < termsForE.length; cti++) {
                    let currentTerm = termsForE[cti];
                    if (currentTerm.length > 0) {
                        var altTerm = [[currentTerm[0], currentTerm[1]]];
                      
                        // check if any of the terms in e's replacement
                        // terms also contains one of the edges to remove
                        if (currentTerm[1] != null) {
                            if (currentTerm[1].some(x => (edgesToRemove.includes(parseInt(x))))) {
                            foundReplacement = false
                      
                            // if so, then we need to replace that term
                            var newTerm = [[1, []]];
                            for (let ttt = 0; ttt < currentTerm[1].length; ttt++) {
                                let tt = currentTerm[1][ttt];
                                if (edgesToRemove.includes(tt)) {
                                    var nt = [];
                                    for (let rdi = 0; rdi < reduceDict[tt].length; rdi++) {
                                        let rd = reduceDict[tt][rdi];
		                        for (let nt1i = 0; nt1i < newTerm.length; nt1i++) {
                                            let nt1 = newTerm[nt1i];
                                            var nt11 = rd[1];
                                            if (nt1[1].length > 0) { 
                                                nt11 = nt1[1].concat(rd[1]);
                                            }
                                            nt.push([parseInt(nt1[0])*parseInt(rd[0]), nt11]);
		            	      }
		            	  }
                                    if (nt.length > 0) { newTerm = nt; }
                                } else {
                                    newTerm = newTerm.map(
                                        function(x) {
                                            if (x[1].length > 0) {
                                                return [x[0], x[1].concat(tt)];
                                            } else {
                                                return [x[0], [tt]];
                                            }
                                    });
		                }
                            }
                            if (newTerm[0][1].length > 0) {
		                altTerm = newTerm.map(x => [x[0]*currentTerm[0], x[1]]);
                            }
                        }}
		        altTermsForE.push(...altTerm);
                    }
		    termsForE = altTermsForE;
                }
            } while (!foundReplacement && (ctr < edgesToRemove.length));
	    reduceDict[e] = termsForE;
        }

        // reduce the potential by replacing each of the terms with its 
        // image in the replacement dictionary.
        var wPrime = [];
        for (let tci = 0; tci < thePotential.length; tci++) {

            let coef = thePotential[tci][0];
            let term = thePotential[tci][1].split(",").map(x => parseInt(x));

            var newTerm = [[coef, []]];
            for (let tti = 0; tti < term.length; tti++) {
                let tt = term[tti];

                var thisLevel = [];
                for (let nt1i = 0; nt1i < newTerm.length; nt1i++) {
                    var nt1 = newTerm[nt1i];

                    for (let rdi = 0; rdi < reduceDict[tt].length; rdi++) {
                        var rd = reduceDict[tt][rdi];
                        thisLevel.push([parseInt(nt1[0])*parseInt(rd[0]),
                                        nt1[1].concat(rd[1])]);
                    }
                }
                if (thisLevel.length > 0) { newTerm = deepCopy(thisLevel); }
            }
            wPrime.push(...newTerm);
        }
        wp = wPrime.map(function(x) {
            return [x[0], x[1].filter(y => y != null)];
        });

	//  make sure all terms are uniquely represented
        // (cycles are written with minimal element first)
	wPrime = [];
	for (let i = 0; i < wp.length; i++) {
            let x = wp[i];
	    let c = x[0];
	    let t = cycleOrder(x[1]).toString();

	    if (wPrime.some(y => y[1] == t)) {
                let idx = wPrime.findIndex(x => (x[1] == t));
                wPrime[idx] = [wPrime[idx][0] + c, t];
            } else if (c != 0) {
                wPrime.push([c, t]);
            }
	}
	// remove all terms with 0 coefficient
	wPrime = wPrime.filter(y => (y[0] != 0));

        return removeEdges(edgesToRemove, QP, altPotential=wPrime);
    } else {
        return deepCopy(QP);
    }
}

function removeEdges(edgeIndices, QP, altPotential="None") {
    var edgesToKeep = [...Array(QP.edges.length).keys()].map(
        function(x) {
            if (edgeIndices.includes(parseInt(x))) {
                return -1;
            } else {
                return parseInt(x);
            }
        });
    var edgeIndexLookup = edgesToKeep.filter(x => x >= 0);
    var edgeIndexLookupBackwards = edgesToKeep.map(function(x) {if (x >= 0) {return edgeIndexLookup.indexOf(x);}});

    // re-index the edges to delete those that are no longer included
    var newEdges = edgeIndexLookup.map(x => QP.edges[x]);

    // update the terms in the potential to use the new edge indexing
    var newPotential = altPotential;
    if (newPotential == "None") {
        newPotential = QP.potential.map(x => [...x]);
    }
    newPotential = newPotential.map(
        function(x){
            let y = x[1].split(",").map(y => edgeIndexLookupBackwards[parseInt(y)]);
	    return [parseInt(x[0]), y.filter(x => x != null).toString()];
        });
    return makeQP(newEdges, QP.nodes, QP.frozenNodes, newPotential, "fromQP");
}
