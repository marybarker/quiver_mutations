var nodes, edges, network;
var frozen_nodes, coefs, terms;

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
        });
    } catch (err) {
        alert(err);
    }
}

function updateNode() {
    try {
        nodes.update({
            id: document.getElementById("node-id").value,
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
        { id: "0"},
        { id: "1"},
        { id: "2"},
        { id: "3"},
        { id: "4"},
        { id: "5"},
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
	  4
      );
  });

  // create a network
  var container = document.getElementById("mynetwork");
  var data = {
      nodes: nodes,
      edges: edges,
  };
  var options = {
      physics: {
          enabled: false,
      },
      interaction: {navigationButtons: true}
  };
  network = new vis.Network(container, data, options);
}
