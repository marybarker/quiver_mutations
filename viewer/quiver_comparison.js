var comparisonContainer = document.getElementById('quiver-comparison-container')

var comparisonExpected = document.getElementById('comparison-expected-input')
var comparisonPotentials = document.getElementById('comparison-potentials-input')

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
  nodes.add(qp.nodes.map(function (x) {
    const i = x.id || x
    return {
      id: i.toString(),
      label: i.toString(),
      x: parseFloat(x.x || (i * 40).toString()),
      y: parseFloat(x.y || (i * 40).toString())
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
  if (!comparisonExpected.value || !comparisonPotentials.value) {
    comparisonContainer.innerHTML = ''
    return
  }

  var expected = JSON.parse(comparisonExpected.value)
  var potentials = JSON.parse(comparisonPotentials.value)

  console.log(expected, potentials)

  var h = document.createElement('h2')
  h.textContent = 'Expected'
  comparisonContainer.appendChild(h)
  var row = document.createElement('div')
  row.className = 'comp-row'
  expected.forEach(function (eQ) {
    var container = document.createElement('div')
    container.className = 'comp-container'
    displayVis(eQ, container)
    row.appendChild(container)
  })
  comparisonContainer.appendChild(row)

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

  potentials.forEach(function (thisPotential, i) {
    // edges and nodes come from the global viz object, potential comes from the text input
    var baseQP = makeQP(edges, nodes, frozen_nodes, potential, 'fromVisDataSet')
    baseQP.potential = thisPotential

    console.log(baseQP)

    var resultQuivers = getAllMutationsForQP(baseQP)

    var h = document.createElement('h2')
    h.textContent = i + ' - ' + JSON.stringify(thisPotential)
    comparisonContainer.appendChild(h)
    var row = document.createElement('div')
    row.className = 'comp-row'
    console.log('found', resultQuivers.quivers.length)
    resultQuivers.quivers.forEach(function (qp) {
      console.log(potential, qp)
      var container = document.createElement('div')
      container.className = 'comp-container'
      displayVis(qp, container)
      row.appendChild(container)
    })
    comparisonContainer.appendChild(row)
  })
}

comparisonExpected.addEventListener('input', performQuiverComparison)
comparisonPotentials.addEventListener('input', performQuiverComparison)
