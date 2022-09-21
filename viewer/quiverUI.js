var nodes, edges, network
var frozen_nodes
var potential
var potentialVisualization
var output_fields = ['id', 'from', 'to', 'coef'] // which data objects to print to screen

document.getElementById('graph-data-edges').addEventListener('change', function () {
  if (this.checked) {
    network.setData({ nodes: nodes, edges: edges })
  }
})

document.getElementById('graph-data-potential').addEventListener('change', function () {
  if (this.checked) {
    network.setData({ nodes: nodes, edges: potentialVisualization })
  }
})

function updateInstructions () {
  const click_mode = document.getElementById('edit-quiver-type').value
  var textOutput = 'Click on canvas to add a node'

  if (click_mode == 'add-edge') {
    textOutput = 'Select a pair of nodes to add an edge between them: control + click to select two nodes'
  } else if (click_mode == 'add-loop') {
    textOutput = 'Click on a node to add a loop'
  } else if (click_mode == 'remove-edge') {
    textOutput = 'Click on an edge to remove it'
  } else if (click_mode == 'add-node') {
    textOutput = 'Click on the canvas to create a node'
  } else {
    textOutput = 'Select a node to perform action'
  }
  document.getElementById('instructions').innerText = textOutput
}

function resolve_click_event (n, p) {
  var click_mode = 'NONE'
  try {
    click_mode = document.getElementById('edit-quiver-type').value
  } catch (err) {
    alert(err)
  }

  if (click_mode == 'add-node') {
    if ((p.nodes.length == 0) && (p.edges.length == 0)) { // make sure we're not superimposing nodes
	    var nv = getUniqueNodeId()
      nodes.add({ id: nv, label: nv, x: p.pointer.canvas.x, y: p.pointer.canvas.y })
    }
  } else if (click_mode == 'remove-node') {
    if (p.nodes.length > 0) {
      var n = p.nodes[0].toString()
	    removeGlobalTerms(nodeIdsToRemove = [n], edgeIdsToRemove = [], potentialIdsToRemove = [])
    }
  } else if (click_mode == 'add-edge') {
    var ne = getUniqueEdgeId()
    if (p.nodes.length > 1) {
	    const n1 = p.nodes[0].toString()
	    const n2 = p.nodes[1].toString()
      edges.add({ id: ne, from: n1, to: n2, arrows: 'to', title: 'edge ' + ne })
    }
  } else if (click_mode == 'remove-edge') {
    if (p.edges.length > 0) {
	    var e = p.edges[0].toString()
	    removeGlobalTerms([], [e], [])
    }
  } else if (click_mode == 'add-loop') {
    var ne = getUniqueEdgeId()
    if (p.nodes.length > 0) {
	    const v = p.nodes[0].toString()
      edges.add({ id: ne, from: v, to: v })
    }
  } else if (click_mode == 'freeze-node') {
    if (p.nodes.length > 0) {
      frozen_nodes.add({ id: p.nodes[0].toString() })
    }
  } else if (click_mode == 'unfreeze-node') {
    if (p.nodes.length > 0) {
      frozen_nodes.remove({ id: p.nodes[0].toString() })
    }
  } else if (click_mode == 'mutate') {
    if (p.nodes.length > 0) {
      const current_node = p.nodes[0].toString()
      var mQP = makeQP(edges, nodes, frozen_nodes, potential)

      mQP = mutateQP(current_node, mQP)
	    updateGlobalsFromQP(mQP)
    }
  }
}

function updateGlobalsFromQP (QP) {
  // NOTE: this routine does not update nodes global.
  var outputs = []
  for (let i = 0; i < QP.edges.length; i++) {
    outputs.push({
	    id: i.toString(),
	    to: QP.edges[i][1].toString(),
	    from: QP.edges[i][0].toString(),
	    arrows: 'to',
	    title: 'edge ' + i.toString()
    })
  }
  edges.clear()
  edges.add(outputs)

  potential.clear()
  for (let i = 0; i < QP.potential.length; i++) {
    addTermToPotential(QP.potential[i][1].toString(), QP.potential[i][0])
  }
}

function updateQPFromJSON (JSONData) {
  clearQP()
  edges.add(JSONData.edges.map(function (x) {
    const i = x.id
    const f = x.from
    const t = x.to
    return {
    	    id: i.toString(),
      from: f.toString(),
      to: t.toString(),
      arrows: 'to',
      title: 'edge ' + i.toString()
    }
  }))
  nodes.add(JSONData.nodes.map(function (x) {
    const i = x.id
    return {
    	    id: i.toString(),
      label: i.toString(),
      x: parseFloat(x.x),
      y: parseFloat(x.y)
    }
  }))
  for (let pi = 0; pi < JSONData.potential.length; pi++) {
    const x = JSONData.potential[pi]
    const i = cycleOrder(x.id.split(',')).toString()
    const c = x.coef
    addTermToPotential(i, c)
  }
  frozen_nodes.add(JSONData.frozenNodes.map(function (x) {
    const i = x.id
    return {
    	id: i.toString()
    }
  }))
}

function removeGlobalTerms (nodeIdsToRemove = [], edgeIdsToRemove = [], potentialIdsToRemove = []) {
  // removes nodes, edges, and potential terms and updates the index/label ordering on edges and nodes

  var oldNodeIds = nodes.getIds()
  var newNodeIndexing = []

  if (nodeIdsToRemove.length > 0) {
    for (let ni = 0; ni < nodeIdsToRemove.length; ni++) {
	    const n = nodeIdsToRemove[ni].toString()
      var extraEdgesToRemove = edges.getIds().filter(x => ((edges.get(x).to == n) || (edges.get(x).from == n)))
      for (let ei = 0; ei < extraEdgesToRemove.length; ei++) {
	        const e = extraEdgesToRemove[ei].toString()
        if (!edgeIdsToRemove.includes(e)) {
          edgeIdsToRemove.push(e)
        }
      }
    }
    // reindex/relabel nodes
    oldNodeIds = nodes.getIds().filter(x => !nodeIdsToRemove.includes(x))
    newNodeIndexing = oldNodeIds.map(function (i) {
      return {
        id: oldNodeIds.indexOf(i).toString(),
        label: oldNodeIds.indexOf(i).toString(),
        x: nodes.get(i).x,
        y: nodes.get(i).y
      }
    })
  }

  var oldEdgeIds = edges.getIds()
  var newEdgeIndexing = []
  var newPotentialTerms = []
  if (edgeIdsToRemove.length > 0) {
    for (let ei = 0; ei < edgeIdsToRemove.length; ei++) {
	    const e = edgeIdsToRemove[ei].toString()
      edges.remove({ id: e })
    }
    // reindex/relabel edges
    oldEdgeIds = edges.getIds()
    newEdgeIndexing = oldEdgeIds.map(function (i) {
	    return {
		    id: oldEdgeIds.indexOf(i).toString(),
	            from: oldNodeIds.indexOf(edges.get(i).from).toString(),
	            to: oldNodeIds.indexOf(edges.get(i).to).toString(),
	            arrows: 'to',
	            title: 'edge ' + oldEdgeIds.indexOf(i).toString()
	        }
	    })

    var potentialTermsToRemove = potential.getIds().filter(x => !potentialTermIsSubsetOfEdges(x))
    // remove terms from potential
    for (let pt = 0; pt < potentialTermsToRemove.length; pt++) {
      potential.remove({ id: potentialTermsToRemove[pt].toString() })
    }
  }
  for (let pt = 0; pt < potentialIdsToRemove.length; pt++) {
    potential.remove({ id: potentialIdsToRemove[pt].toString() })
  }

  // relabel edges in potential by their new numbering
  for (let pt = 0; pt < potential.getIds().length; pt++) {
    const t = potential.getIds()[pt]
    newPotentialTerms.push({ id: t.split(',').map(x => oldEdgeIds.indexOf(x.toString())).toString(), coef: potential.get(t).coef })
  }
  var newFn = frozen_nodes.getIds().filter(x => oldNodeIds.includes(x)).map(x => { id: oldNodeIds.indexOf(x).toString() })
  if (newNodeIndexing.length > 0) {
    nodes.clear()
    nodes.add(newNodeIndexing)
  }
  if (newEdgeIndexing.length > 0) {
    edges.clear()
    edges.add(newEdgeIndexing)
  }
  if (newPotentialTerms.length > 0) {
    potential.clear()
    potential.add(newPotentialTerms)
  }
  if (newFn.legnth > 0) {
    frozen_nodes.clear()
    frozen_nodes.add(newFn)
  }
}

function potentialTermIsSubsetOfEdges (term) {
  var esInTerm = term.split(',')
  var currentEdges = edges.getIds()
  return (esInTerm.filter(x => !currentEdges.includes(x.toString())).length < 1)
}

function addTermToPotential (t, coef = 1) {
  var c2 = 0
  try {
    c2 = potential.get(t)
  } catch (err) {}
  if (potentialTermIsSubsetOfEdges(t)) {
    if (c2 == null) {
      potential.add({ id: t, coef: coef.toString() })
    } else {
      const c3 = parseFloat(coef) + parseFloat(c2.coef)
      if (c3 > 0 || c3 < 0) {
        potential.update({ id: t, coef: c3.toString() })
      } else {
        potential.remove({ id: t })
      }
    }
  } else {
    alert('Error updating potential: term ' + t + ' contains an invalid edge')
  }
}

function generateRandomPotential () {
  var np = randomPotential(nodes, edges)

  potential.clear()
  for (let i = 0; i < np.length; i++) {
    const x = np[i]
    if (x[1] != ',') {
      addTermToPotential(x[1], x[0])
    }
  }
  /*
      potential.add(np.filter(x => (x[1] != ",")).map(function(x) {
      return {
          id: x[1],
          coef: x[0].toString(),
          };
      }));
      */
}

function getUniqueEdgeId () { /* create a string value that is not currently in ids for edges */
  var ne = edges.getIds()
  ne = ne.map(function (x) { return parseInt(x) })
  return ne.reduce(function (a, b) { return Math.max(a, b) + 1 }, 0).toString()
}

function getUniqueNodeId () { /* create a string value that is not currently in ids for nodes */
  var nv = nodes.getIds()
  nv = nv.map(function (x) { return parseInt(x) })
  return nv.reduce(function (a, b) { return Math.max(a, b) + 1 }, 0).toString()
}

function clearPotential () {
  potential.clear()
}

function clearQP () {
  edges.clear()
  frozen_nodes.clear()
  nodes.clear()
  potential.clear()
}

function updatePotential () {
  var t = document.getElementById('term-input').value
  var c1 = parseFloat(document.getElementById('coefficient-input').value)
  addTermToPotential(t, c1)
}

// assigns each self-loop a unique radius so that they don't overlap
function updateEdgeRadii (edges, id) {
  var thisEdge = edges.get(id)

  if (thisEdge.from === thisEdge.to) {
    var count = edges.get().filter(function (otherEdge) {
      return otherEdge.from === thisEdge.from & otherEdge.to === thisEdge.to && parseInt(otherEdge.id) < parseInt(thisEdge.id)
    }).length

    thisEdge.selfReference = {
      size: 15 + (count * 5)
    }

    edges.update(thisEdge)
  }
}

function draw () {
  // create an array with nodes
  nodes = new vis.DataSet()
  nodes.on('*', function () {
    document.getElementById('nodes').innerText = JSON.stringify(
      nodes.get(),
      output_fields,
      4
    )
  })
  nodes.add([
    { id: '0', label: '0', x: -100, y: 0 },
    { id: '1', label: '1', x: 0, y: 100 },
    { id: '2', label: '2', x: 100, y: 0 },
    { id: '3', label: '3', x: 100, y: -100 },
    { id: '4', label: '4', x: 0, y: -200 },
    { id: '5', label: '5', x: -100, y: -100 }
  ])

  // create an array with edges
  edges = new vis.DataSet()
  edges.on('*', function () {
    document.getElementById('edges').innerText = JSON.stringify(
      edges.get(),
      output_fields,
      4
    )
  })
  edges.add([
    { id: '0', from: '0', to: '1', arrows: 'to', title: 'edge 0' },
    { id: '1', from: '1', to: '0', arrows: 'to', title: 'edge 1' },
    { id: '2', from: '1', to: '2', arrows: 'to', title: 'edge 2' },
    { id: '3', from: '2', to: '1', arrows: 'to', title: 'edge 3' },
    { id: '4', from: '0', to: '2', arrows: 'to', title: 'edge 4' },
    { id: '5', from: '2', to: '0', arrows: 'to', title: 'edge 5' },
    { id: '6', from: '2', to: '3', arrows: 'to', title: 'edge 6' },
    { id: '7', from: '3', to: '2', arrows: 'to', title: 'edge 7' },
    { id: '8', from: '3', to: '3', arrows: 'to', title: 'edge 8' },
    { id: '9', from: '3', to: '4', arrows: 'to', title: 'edge 9' },
    { id: '10', from: '4', to: '3', arrows: 'to', title: 'edge 10' },
    { id: '11', from: '4', to: '4', arrows: 'to', title: 'edge 11' },
    { id: '12', from: '4', to: '5', arrows: 'to', title: 'edge 12' },
    { id: '13', from: '5', to: '4', arrows: 'to', title: 'edge 13' },
    { id: '14', from: '5', to: '5', arrows: 'to', title: 'edge 14' },
    { id: '15', from: '5', to: '5', arrows: 'to', title: 'edge 15' },
    { id: '16', from: '5', to: '0', arrows: 'to', title: 'edge 16' },
    { id: '17', from: '0', to: '5', arrows: 'to', title: 'edge 17' }
  ])

  // update the initial dataset

  edges.get().forEach(edge => updateEdgeRadii(edges, edge.id))

  // and whenever an edge is added
  edges.on('add', function (event, properties, senderId) {
    properties.items.forEach(function (i) {
      updateEdgeRadii(edges, i)
    })
  })

  frozen_nodes = new vis.DataSet()
  frozen_nodes.on('*', function () {
    document.getElementById('frozen_nodes').innerText = JSON.stringify(
      frozen_nodes.get(),
      output_fields,
      4)
  })
  potential = new vis.DataSet()
  potential.on('*', function () {
    document.getElementById('potential').innerText = JSON.stringify(
      potential.get(),
      output_fields,
      4)
  })

  potentialVisualization = new vis.DataSet()
  potential.on('*', function () {
    potentialVisualization.clear()
    potential.get().forEach(function (term, idx) {
      var termEdges = term.id.split(',')
      var termColor = 'rgb(' + [Math.round(50 + Math.random() * 150), Math.round(50 + Math.random() * 150), Math.round(50 + Math.random() * 150)].join(',') + ')'
      for (var i = 0; i < termEdges.length; i++) {
        potentialVisualization.add(Object.assign({}, edges.get(termEdges[i]), {
          id: null, // need to remove the original edge ID so that each item here gets assigned a new unique ID
          color: termColor,
          label: term.coef,
          font: { size: 40 },
          selfReference: { size: 40 + (5 * idx) }
        }))
      }
    })
  })

  potential.add([
    { id: '0,2,5', coef: '1' },
    { id: '1,4,3', coef: '1' },
    { id: '8,9,10', coef: '1' },
    { id: '2,6,7,3', coef: '1' },
    { id: '0,1,17,16', coef: '1' },
    { id: '6,8,7', coef: '1' }
  ])

  // create a network
  var container = document.getElementById('mynetwork')
  var data = {
    nodes: nodes,
    edges: edges
  }
  var options = {
    interaction: { hover: true, multiselect: true },
    nodes: {
      borderWidth: 1,
      size: 45,
      color: {
        border: '#222222',
        background: 'grey'
      },
      font: {
        color: 'black',
        size: 11,
        face: 'arial'
      },
      physics: false
    },
    edges: {
      arrows: {
        to: { enabled: true }
      },
      color: {
        color: '#848484',
        highlight: '#848484',
        hover: '#848484'
      },
      font: {
        color: '#343434',
        size: 11, // px
        face: 'arial',
        background: 'none',
        strokeWidth: 5, // px
        strokeColor: '#ffffff'
      }
      // physics: {enabled:true},
    }
  }
  network = new vis.Network(container, data, options)
  // network.setOptions({physics:{enabled:false}});
  network.on('click', function (params) {
    resolve_click_event(network, params)
  })
}

function showExchangeNumber () {
  const output = document.getElementById('exchange-number-output')
  try {
    const result = getAllMutationsForQP(makeQP(edges, nodes, frozen_nodes, potential))
    if (result.timeout) {
      output.textContent = 'Timed out'
    } else {
      output.textContent = result.quivers.length
    }
  } catch (e) {
    console.error(e)
    output.textContent = 'Error'
  }
}
