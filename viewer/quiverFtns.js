var nodes, edges, network;
var frozen_nodes;
var potential;
var editMode = "Quiver"

function resolve_click_event(n, p) {
    //if (editMode == "Quiver") {
	// first identify what type of editing option is enabled(adding/removing/freezing)
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
        }
    //}
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
function updateEditMode(editModeOption) {
    editMode = editModeOption.value;
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
function toJSON(obj) {
    return JSON.stringify(obj, null, 4);
}
function freezeNode() {
    try {
        frozen_nodes.add({ id: document.getElementById("node-id").value });
    } catch (err) {
        alert(err)
    }
}
function unfreezeNode() {
    try {
        frozen_nodes.remove({ id: document.getElementById("node-id").value });
    } catch (err) {
        alert(err)
    }
}
function addNode() {
    try {
        nodes.add({
            id: document.getElementById("node-id").value,
            label: document.getElementById("node-id").value,
        });
    } catch (err) {
        alert(err);
    }
}
function updateNode() {
    try {
        nodes.update({
            id: document.getElementById("node-id").value,
            label: document.getElementById("node-id").value,
        });
    } catch (err) {
        alert(err);
    }
}
function removeNode() {
    try {
        nodes.remove({ id: document.getElementById("node-id").value });
    } catch (err) {
        alert(err);
    }
}
function addEdge() {
    try {
        edges.add({
            id: document.getElementById("edge-id").value,
            from: document.getElementById("edge-from").value,
            to: document.getElementById("edge-to").value,
            arrows: "to",
        });
    } catch (err) {
        alert(err);
    }
}
function updateEdge() {
    try {
        edges.update({
            id: document.getElementById("edge-id").value,
            from: document.getElementById("edge-from").value,
            to: document.getElementById("edge-to").value,
            arrows: "to",
        });
    } catch (err) {
        alert(err);
    }
}
function removeEdge() {
    try {
        edges.remove({ id: document.getElementById("edge-id").value });
    } catch (err) {
        alert(err);
    }
}
function draw() {
    // create an array with nodes
    nodes = new vis.DataSet();
    nodes.on("*", function () {
        document.getElementById("nodes").innerText = JSON.stringify(
            nodes.get(),
            null,
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
            null,
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
          null,
	  4);
    });
    potential = new vis.DataSet();
    potential.on("*", function () {
    document.getElementById("potential").innerText = JSON.stringify(
          potential.get(),
          null,
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
