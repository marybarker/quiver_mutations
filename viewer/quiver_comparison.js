var comparisonContainer = document.getElementById('quiver-comparison-container')

var comparisonExpected = document.getElementById('comparison-expected-input')
var comparisonPotentials = document.getElementById('comparison-potentials-input')

/*
The array of all possible quivers has inconsistent node IDs
This uses the position info to remap them so the node IDs all match those of the first quiver
*/
function remapQPNodes (allQPSArr) {
  return deepCopy(allQPSArr).map(function (qp, idx) {
    if (idx !== 0) {
      var oldNewMap = {}
      qp.nodes.forEach(function (node) {
        const newId = allQPSArr[0].nodes.find(n => Math.abs(n.x - node.x) < 0.000001 && Math.abs(n.y - node.y) < 0.000001).id
        oldNewMap[node.id] = newId
        node.id = newId
      })
      qp.edges.forEach(function (edge, idx) {
        edge.from = oldNewMap[edge.from]
        edge.to = oldNewMap[edge.to]
      })
    }
    return qp
  })
}

function displayVis (qp, container) {
  var nodes = new vis.DataSet()
  var edges = new vis.DataSet()
  var frozen_nodes = new vis.DataSet()

  // from updateqpfromjson
  edges.add(qp.edges.map(function (x, idx) {
    const i = x.id || idx.toString()
    const f = x.from || x[0].toString()
    const t = x.to || x[1].toString()
    return {
      id: i.toString(),
      from: f.toString(),
      to: t.toString(),
      arrows: 'to',
      title: 'edge ' + i.toString()
    }
  }))

  edges.get().forEach(edge => updateEdgeRadii(edges, edge.id))

  // and whenever an edge is added
  edges.on('add', function (event, properties, senderId) {
    properties.items.forEach(function (i) {
      updateEdgeRadii(edges, i)
    })
  })

  nodes.add(qp.nodes.map(function (x, i) {
    const id = x.id || x
    return {
      id: id.toString(),
      label: id.toString(),
      x: (x.x || 25 * i).toString(),
      y: (x.y || 25 * i).toString()
    }
  }))

  frozen_nodes.add(qp.frozenNodes.map(function (x) {
    const i = x.id || x
    return {
      id: i.toString()
    }
  }))

  var data = {
    nodes: nodes,
    edges: edges
  }
  var options = {
    clickToUse: true,
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
  var network = new vis.Network(container, data, options)
}

function performQuiverComparison () {
  // string representing the edge counts at each node
  function stringifyQuiver (quiver) {
    // need to count loops and regular edges separately
    // otherwise, two nodes A and B, both with self loops are counted the same as having two edges connecting A and B
    var edgeCounts = new Array(quiver.nodes.length)
    var loopCounts = new Array(quiver.nodes.length).fill(0)

    // can't use array.fill here because it doesn't create unique array instances
    for (var i = 0; i < edgeCounts.length; i++) {
      edgeCounts[i] = [0, 0]
    }

    quiver.edges.forEach(function (e) {
      if (Object.hasOwn(e, 'from')) {
        if (e.from === e.to) {
          loopCounts[e.from]++
        } else {
          edgeCounts[parseInt(e.from)][0]++
          edgeCounts[parseInt(e.to)][1]++
        }
      } else {
        if (e[0] === e[1]) {
          loopCounts[e[0]]++
        } else {
          edgeCounts[parseInt(e[0])][0]++
          edgeCounts[parseInt(e[1])][1]++
        }
      }
    })

    var stringEdges = edgeCounts.map(v => v[0] + '.' + v[1]).concat(loopCounts)

    return stringEdges.sort().join(',')
  }

  comparisonContainer.innerHTML = ''

  if (!comparisonExpected.value || !comparisonPotentials.value) {
    return
  }

  var expected = remapQPNodes(JSON.parse(comparisonExpected.value))
  var potentials = JSON.parse(comparisonPotentials.value)

  var toprow = document.createElement('div')
  toprow.id = 'comparison-top-row'
  var h = document.createElement('h3')
  h.textContent = 'Expected'
  toprow.appendChild(h)
  var row = document.createElement('div')
  row.className = 'comp-row'
  expected.forEach(function (eQ) {
    var container = document.createElement('div')
    container.className = 'comp-container'
    displayVis(eQ, container)
    row.appendChild(container)
  })
  toprow.appendChild(row)
  comparisonContainer.appendChild(toprow)

  var expectedStrings = expected.map(q => stringifyQuiver(q))

  potentials.sort((a, b) => a.length - b.length).slice(0, 20).forEach(function (thisPotential, i) {
    const base = convertQuiver(deepCopy(expected[0]))
    var baseQP = makeQP(base.edges, base.nodes, base.frozenNodes, base.potential, 'fromThing')

    baseQP.potential = thisPotential

    var resultQuivers = getAllMutationsForQP(baseQP).quivers.map(q => JSON.parse(q))

    resultQuivers = resultQuivers.map(function (resQ) {
      const expectedBase = deepCopy(expected[0])
      expectedBase.edges = resQ.edges.map(function (edg, i) {
        return {
          id: i.toString(),
          to: edg[1].toString(),
          from: edg[0].toString(),
          arrows: 'to',
          title: 'edge ' + i.toString()
        }
      })
      return expectedBase
    })

    resultQuivers = resultQuivers.sort((a, b) => {
      return expectedStrings.indexOf(stringifyQuiver(a)) - expectedStrings.indexOf(stringifyQuiver(b))
    })

    // is this an exact match?

    const expectedMatchStrings = deepCopy(expected).map(q => stringifyQP(convertQuiver(q)))
    const actualMatchStrings = deepCopy(resultQuivers).map(q => stringifyQP(convertQuiver(q)))

    const exactMatch = arrayEquals(expectedMatchStrings.sort(), actualMatchStrings.sort())

    console.log(expectedMatchStrings.sort(), actualMatchStrings.sort())

    var h = document.createElement('h3')
    h.textContent = (exactMatch ? '[exact match] ' : '') + i + ' - ' + JSON.stringify(thisPotential)
    comparisonContainer.appendChild(h)
    var row = document.createElement('div')
    row.className = 'comp-row'
    resultQuivers.forEach(function (qp) {
      var container = document.createElement('div')
      container.className = 'comp-container'
      displayVis(qp, container)
      row.appendChild(container)
    })
    comparisonContainer.appendChild(row)
  })

  if (potentials.length > 20) {
    var h = document.createElement('h1')
    h.textContent = (potentials.length - 20) + ' more hidden'
    comparisonContainer.appendChild(h)
  }
}

comparisonExpected.addEventListener('input', performQuiverComparison)
comparisonPotentials.addEventListener('input', performQuiverComparison)
