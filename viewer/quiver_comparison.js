var comparisonContainer = document.getElementById('quiver-comparison-container')

var comparisonExpected = document.getElementById('comparison-expected-input')
var comparisonPotentials = document.getElementById('comparison-potentials-input')

function displayVis (qp, container) {
  var nodes = new vis.DataSet()
  var edges = new vis.DataSet()
  var frozen_nodes = new vis.DataSet()

  function circleCoord (i, n) {
    return {
      x: 125 + (125 * Math.cos((2 * Math.PI) * (i / n))),
      y: 125 + (125 * Math.sin((2 * Math.PI) * (i / n)))
    }
  }

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

  nodes.add(qp.nodes.map(function (x) {
    const i = x.id || x
    return {
      id: i.toString(),
      label: i.toString(),
      x: circleCoord(i, qp.nodes.length).x.toString(),
      y: circleCoord(i, qp.nodes.length).y.toString()
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
    interaction: { hover: true },
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
 	   physics: { enabled: false }
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
        strokeColor: '#ffffff',
        align: 'vertical'
      }
 	   // physics: {enabled:true},
    },
    interaction: { multiselect: true },
    navigation: true
  }
  var network = new vis.Network(container, data, options)
}

function performQuiverComparison () {
  // string representing the edge counts at each node
  function stringifyQuiver (quiver) {
    var edgeCounts = new Array(quiver.nodes.length)
    // can't use array.fill here because it doesn't create unique array instances
    for (var i = 0; i < edgeCounts.length; i++) {
      edgeCounts[i] = [0, 0]
    }

    quiver.edges.forEach(function (e) {
      if (Object.hasOwn(e, 'from')) {
        edgeCounts[parseInt(e.from)][0]++
        edgeCounts[parseInt(e.to)][1]++
      } else {
        edgeCounts[parseInt(e[0])][0]++
        edgeCounts[parseInt(e[1])][1]++
      }
    })

    var stringEdges = edgeCounts.map(v => v[0] + '.' + v[1])

    return stringEdges.sort().join(',')
  }

  comparisonContainer.innerHTML = ''

  if (!comparisonExpected.value || !comparisonPotentials.value) {
    return
  }

  var expected = JSON.parse(comparisonExpected.value)
  var potentials = JSON.parse(comparisonPotentials.value)

  console.log(expected, potentials)

  var toprow = document.createElement('div')
  toprow.id = 'comparison-top-row'
  var h = document.createElement('h2')
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

  function convertQuiver (q) {
    for (var i = 0; i < q.edges.length; i++) {
      if (!(q.edges[i] instanceof Array)) {
        q.edges[i] = [parseInt(q.edges[i].from), parseInt(q.edges[i].to)]
      }
    }
    for (var i = 0; i < q.nodes.length; i++) {
      if (q.nodes[i] instanceof Object) {
        q.nodes[i] = parseInt(q.nodes[i].id)
      }
    }

    for (var i = 0; i < q.frozenNodes.length; i++) {
      if (q.frozenNodes[i] instanceof Object) {
        q.frozenNodes[i] = parseInt(q.frozenNodes[i].id)
      }
    }

    return q
  }

  var expectedStrings = expected.map(q => stringifyQuiver(q))

  potentials.sort((a, b) => a.length - b.length).slice(0, 20).forEach(function (thisPotential, i) {
    // edges and nodes come from the global viz object, potential comes from the text input
    var baseQP = makeQP(edges, nodes, frozen_nodes, potential, 'fromVisDataSet')
    baseQP.potential = thisPotential

    console.log(baseQP)

    var resultQuivers = getAllMutationsForQP(baseQP)

    resultQuivers.quivers = resultQuivers.quivers.sort((a, b) => {
      return expectedStrings.indexOf(stringifyQuiver(a)) - expectedStrings.indexOf(stringifyQuiver(b))
    })

    console.log(expectedStrings, resultQuivers.quivers.map(q => stringifyQuiver(q)))

    var h = document.createElement('h2')
    h.textContent = i + ' - ' + JSON.stringify(thisPotential)
    comparisonContainer.appendChild(h)
    var row = document.createElement('div')
    row.className = 'comp-row'
    console.log('found', resultQuivers.quivers.length)
    resultQuivers.quivers.forEach(function (qp) {
      console.log(thisPotential, qp)
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
