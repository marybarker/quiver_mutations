var nodes, edges, network;
var frozen_nodes;
var potential;
var editMode = "Quiver"
var output_fields = ["id", "from", "to", "coef"] // which data objects to print to screen

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
            for (let i = 0; i < mQP.edges.length; i++) {
                let o = {id: i.toString(), to: mQP.edges[i][1].toString(), from: mQP.edges[i][0].toString(), arrows: "to"};
                if (i.toString() in edges.getIds()) {
                    edges.update(o);
                } else {
                    edges.add(o);
                }
            }
	    //potential.update(mQP.ps);
        }
    }
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
        potential.add({id: t, coef: c1});
    } else {
        potential.update({id: t, coef: c1+c2.coef});
    }
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

// QP FUNCTIONS
function deepCopy(A) {
    return JSON.parse(JSON.stringify(A));
}


function makeQP(es, ns, fn, p, inputType="fromVisDataSet") {
    // create graph info (arrow to node structures)
    var awh = Array.from(Array(ns.length), x => []);
    var awt = Array.from(Array(ns.length), x => []);
    var la = Array.from(Array(ns.length), x => []);
    var fns = [];
    var theseNodes = [];
    var theseEdges = [];

    if (inputType == "fromVisDataSet") {
        fns = fn.getIds().map(x => parseInt(x));
        theseNodes = ns.getIds().map(x => parseInt(x));
        theseEdges = es.getIds().map(x => [parseInt(es.get(x).from), parseInt(es.get(x).to)]);
    } else {
        fns = Array.from(fns, x => parseInt(x));
        theseNodes = Array.from(ns, x => parseInt(x));
        theseEdges = deepCopy(es).map(x => [parseInt(x[0]), parseInt(x[1])]);
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
        edges: deepCopy(theseEdges),
        frozenNodes: deepCopy(fns),
        loopsAt: la,
        nodes: deepCopy(theseNodes),
        potential: deepCopy(p)
    }
}


function mutateQP(vertex, QP) {
    const v = parseInt(vertex);
    if (QP.canMutate[v]) {
        // reverse the arrows incident to vertex
        var savedEdges = deepCopy(QP.edges);
        for (let i in QP.arrowsWithHead[v]) {
            let e = savedEdges[i];
            savedEdges[parseInt(i)] = [e[1],e[0]];
        }
        for (let i in QP.arrowsWithTail[v]) {
            let e = savedEdges[i];
            savedEdges[parseInt(i)] = [e[1],e[0]];
        }

        var wPrime = [];
        var edgeCtr = QP.edges.length;
        var shortcuts = [];
        // now add 'shortcuts'
        for (let ei1 in QP.arrowsWithHead[v]) {
            var i = QP.edges[parseInt(ei1)][0];
            for (let ei2 in QP.arrowsWithTail[v]) {
                var j = QP.edges[ei2][1];

                shortcuts.push([ei1,ei2]);
                savedEdges.push([i,j]);
                wPrime.push([1, ei2.toString()+","+ei1.toString()+","+edgeCtr.toString()]);
                edgeCtr++;
            }
        }

        // update the potential
        for (let mc in QP.potential) {
            var coef = mc[0];
            var monoid = mc[1];

            var m = " ";
            var foundMatch = false;
            for (let i = 0; i < monoid.length; i++) {
               m1 = parseInt(monoid[i]);
               m2 = parseInt(monoid[(i+1)%monoid.length]);

               const isIn = shortcuts.findIndex(e => e == [m1,m2]);
               if((isIn > 0) && !foundMatch) {
                   var val = QP.edges.length + parseInt(isIn);
                   m = m + val.toString();
                   foundMatch = true;
               } else {
                   m = m + m1.toString();
                   foundMatch = false;
               }
            }
            wPrime.push([coef, m]);
        }
        // reduce the resulting quiver
        //return reduce(makeQP(savedEdges, QP.nodes, QP.frozenNodes, wPrime, inputType="fromQP"));
        var retVal = makeQP(savedEdges, QP.nodes, QP.frozenNodes, wPrime, inputType="fromQP");
        return retVal;
    } else {
        return makeQP(QP.edges, QP.nodes, QP.frozenNodes, QP.potential, inputType="fromQP");
    }
}
