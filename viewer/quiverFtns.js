if (typeof module !== 'undefined') {
  var {unique} = require('./sharedFtns.js')
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*        Begin QP manipulation functions (no global variables)          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
function allThreeCycles (QP) {
  var triples = []
  const alreadySeen = []
  for (const v1 in QP.nodes) {
    const arrowsOut1 = QP.arrowsWithTail[v1]

    for (let e1ii = 0; e1ii < arrowsOut1.length; e1ii++) {
      const e1i = arrowsOut1[e1ii]
      const e1 = QP.edges[e1i]
      const v2 = e1[1]
      const arrowsOut2 = QP.arrowsWithTail[v2]

      for (let e2ii = 0; e2ii < arrowsOut2.length; e2ii++) {
        const e2i = arrowsOut2[e2ii]

        if (e2i != e1i) {
          const e2 = QP.edges[e2i]
          const v3 = e2[1]
          const arrowsOut3 = QP.arrowsWithTail[v3]

          for (let e3ii = 0; e3ii < arrowsOut3.length; e3ii++) {
            const e3i = arrowsOut3[e3ii]

            if ((e3i != e1i) && (e3i != e2i)) {
              const e3 = QP.edges[e3i]
              const v4 = e3[1]

              if (v4 == v1) {
                const triple = cycleOrder([e1i, e2i, e3i]).toString()
                if (!alreadySeen.includes(triple)) {
                  triples.push([e1i, e2i, e3i])
                  alreadySeen.push([e1i, e2i, e3i].toString())
                }
              }
            }
          }
        }
      }
    }
  }
  return triples
}

function arrayEquals (a, b) {
  if (a.length !== b.length) {
    return false
  }
  for (var i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) {
      return false
    }
  }
  return true
}

function combineLikeTermsInPotential (potential) {
  toRet = []
  addedTerms = []
  for (let i = 0; i < potential.length; i++) {
    const termcoef = potential[i]
    const trm = cycleOrder(termcoef[1].split(',')).toString()
    const idx = toRet.map(function (tc, i) { if (tc[1] == trm) { return i } }).filter(x => x != null)[0]
    if (idx != null) {
      toRet[idx] = [toRet[idx][0] + termcoef[0], trm]
    } else {
	    toRet.push([Number(termcoef[0]), trm])
    }
  }
  return toRet.filter(y => Math.abs(y[0]) > 0)
}

function cycleOrder (cycle) {
  // order a list of integers with minimal element first
  // (note: need to fix this for multiple instances of min element)
  const thisCycle = cycle.filter(y => y != null).map(x => parseInt(x))
  const minVal = Math.min(...thisCycle)
  const minIdx = thisCycle.indexOf(minVal)
  return thisCycle.slice(minIdx).concat(thisCycle.slice(0, minIdx))
}

function deepCopy (A) {
  return JSON.parse(JSON.stringify(A))
}

/* clones nested arrays that only contain simple objects (strings and numbers)
A bit faster than deepCopy
https://stackoverflow.com/a/37503916
 */
function simpleArrayClone (arr) {
  return arr.map(e => Array.isArray(e) ? simpleArrayClone(e) : e);
};

function edgesOutOf (vertex, edge_list) {
  return Array.from(Array(edge_list.length).keys()).map(x => parseInt(x)).filter(x => edge_list[x][0] == vertex)
}

function findCycleDFS (begin_vertex, at_vertex, edge_list, edges_so_far, min_length = 0) {
  const edgesOut = edgesOutOf(at_vertex, edge_list).filter(y => !edges_so_far.includes(y))

  for (let ei = 0; ei < edgesOut.length; ei++) {
    const e = edge_list[edgesOut[ei]]
    const esMet = deepCopy(edges_so_far).map(y => parseInt(y))
    esMet.push(edgesOut[ei])
    const currentAt = e[1]

    if ((begin_vertex == currentAt) && (esMet.length >= min_length)) {
      return [true, esMet]
    }
    return findCycleDFS(begin_vertex, currentAt, edge_list, esMet, min_length)
  }
  return [false, esMet]
}

function findDependencies (currentItem, met, lookups) {
  /* find the items that currentItem depends on in the lookup table
     * (this is recursively done, so the tree of dependencies is listed exhaustively)
     * stopping criterion: no new dependencies found.
     * circularity is dealt with by only adding new dependants */
  if (met.includes(currentItem)) {
    return [true, met]
  } else {
    const newMet = met.concat(currentItem)
    const currentLookup = lookups[parseInt(currentItem)]
    for (let i = 0; i < currentLookup.length; i++) {
	    const item = currentLookup[i]
	    return findDependencies(item, newMet, lookups)
    }
    return [false, newMet]
  }
}

function makeQP (es, ns, fn, p, inputType = 'fromVisDataSet') {
  // create graph info (arrow to node structures)
  var awh = Array.from(Array(ns.length), x => [])
  var awt = Array.from(Array(ns.length), x => [])
  var la = Array.from(Array(ns.length), x => [])
  var fns = []
  var theseNodes = []
  var theseEdges = []
  var thisPotential = []

  if (inputType == 'fromVisDataSet') {
    fns = fn.getIds().map(x => parseInt(x))
    theseNodes = ns.getIds().map(x => parseInt(x))
    theseEdges = es.getIds().map(x => [parseInt(es.get(x).from), parseInt(es.get(x).to)])
    // make sure that the edges in each cycle of the potential are now listed by their index,
    // rather than id, since adding/subtracting edges can leave gaps in the id#s
    var edgeIDMap = es.getIds()
    thisPotential = p.getIds().map(x => [parseFloat(p.get(x).coef), x.split(',').map(y => edgeIDMap.indexOf(y)).toString()])
  } else {
    fns = Array.from(fn, x => parseInt(x))
    theseNodes = Array.from(ns, x => parseInt(x))
    theseEdges = es.filter(x => (x != null)).map(x => [parseInt(x[0]), parseInt(x[1])])
    thisPotential = p
  }

  for (let ei = 0; ei < theseEdges.length; ei++) {
    const e = theseEdges[ei]
    if (e[0] == e[1]) {
      la[theseNodes.indexOf(e[0])].push(ei)
    } else {
      awh[theseNodes.indexOf(e[1])].push(ei)
      awt[theseNodes.indexOf(e[0])].push(ei)
    }
  }
  // can_mutate list
  var cm = theseNodes.map(
    function (v) {
      return ((la[v].length < 1) && !(v in fns))
    })
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

function stringifyQP (qp, includePotential = false) {
  /*
  AWT/AWH contain the same information that is already included in edges, so it doesn't need to be included
  And the sorting of edges makes the edge IDs contained in AWH/AWT incorrect, so they need to be removed
  */
  const qpCopy = {
    edges: deepCopy(qp.edges),
    nodes: deepCopy(qp.nodes),
    frozenNodes: deepCopy(qp.frozenNodes),
  }

  if (includePotential) {
    qpCopy.potential = deepCopy(qp.potential)
  }

  qpCopy.nodes.sort()
  qpCopy.frozenNodes.sort()

  qpCopy.edges.sort(function(a, b) {
    if (a[0] === b[0]) {
      return a[1] - b[1]
    } else {
      return a[0] - b[0]
    }
  })

  if (qpCopy.potential) {
    qpCopy.potential.sort(function (a, b) {
      var as = JSON.stringify(a)
      var bs = JSON.stringify(b)
      if (as < bs) {
        return -1
      }
      if (as > bs) {
        return 1
      }
      return 0
    })
  }

  return JSON.stringify(qpCopy)
}

function mutateQP (vertex, QP) {
  const v = parseInt(vertex)
  if (QP.canMutate[v]) {
    // reverse the arrows incident to vertex
    var savedEdges = deepCopy(QP.edges)
    for (const ii in QP.arrowsWithHead[v]) {
      const i = QP.arrowsWithHead[v][ii]
      const e = savedEdges[i]
      savedEdges[parseInt(i)] = [e[1], e[0]]
    }
    for (const ii in QP.arrowsWithTail[v]) {
      const i = QP.arrowsWithTail[v][ii]
      const e = savedEdges[i]
      savedEdges[parseInt(i)] = [e[1], e[0]]
    }

    var delta = []
    var edgeCtr = QP.edges.length
    var shortcuts = []
    // now add 'shortcuts'
    for (const ei1i in QP.arrowsWithHead[v]) {
      const ei1 = QP.arrowsWithHead[v][ei1i]
      var i = QP.edges[parseInt(ei1)][0]
      for (const ei2i in QP.arrowsWithTail[v]) {
        const ei2 = QP.arrowsWithTail[v][ei2i]
        var j = QP.edges[ei2][1]

        shortcuts.push([ei1, ei2])
        savedEdges.push([i, j])
        delta.push([1, ei2.toString() + ',' + ei1.toString() + ',' + edgeCtr.toString()])
        edgeCtr++
      }
    }

    var wPrime = []
    // update the potential
    for (let mci = 0; mci < QP.potential.length; mci++) {
      const mc = QP.potential[mci]
      var coef = mc[0]
      var monoid = mc[1].split(',')
      const ml = monoid.length
      var m = ''
      var foundMatch = false
      for (let i = 0; i < monoid.length; i++) {
        m1 = parseInt(monoid[i])
        m2 = parseInt(monoid[(i + 1) % ml])
        m0 = parseInt(monoid[(i + ml - 1) % ml])

        const isIn = shortcuts.findIndex(e => arrayEquals(e, [m1, m2]))
        const wasIn = shortcuts.findIndex(e => arrayEquals(e, [m0, m1]))
        if ((i > 0) || (wasIn < 0)) {
          if ((isIn >= 0) && !foundMatch) {
            var val = QP.edges.length + parseInt(isIn)
            m = m + val.toString() + ','
            foundMatch = true
          } else {
            if (!foundMatch) {
              m = m + m1.toString() + ','
            }
            foundMatch = false
          }
        }
      }
      wPrime.push([coef, m.slice(0, -1)])
    }
    wPrime.push(...delta)

    // reduce the resulting quiver
    return reduce(makeQP(savedEdges, QP.nodes, QP.frozenNodes, wPrime, inputType = 'fromQP'))
  } else {
    return makeQP(QP.edges, QP.nodes, QP.frozenNodes, QP.potential, inputType = 'fromQP')
  }
}

function getAllMutationsForQP (qp, maxMutationsToFind = Infinity) {
  var alreadySeen = [stringifyQP(qp)]
  var chains = ['']

  var maxRuntime = 10000
  var beginTime = Date.now()

  function collectMutations (qp, chain) {
    for (var i = 0; i < qp.nodes.length; i++) {
      if (!qp.canMutate[i]) {
        continue
      }
      const currentNode = qp.nodes[i]
      if (Date.now() - beginTime > maxRuntime) {
        return
      }
      if (chains.length === maxMutationsToFind) {
        return
      }
      var mutated = mutateQP(currentNode, qp)
      //we don't need a deep copy here, because we don't use qp again after this point (and mutateQP shouldn't modify it anyway)
      var mutatedStr = stringifyQP(mutated)
      if (!alreadySeen.includes(mutatedStr)) {
        alreadySeen.push(mutatedStr)
        chains.push(chain + currentNode)
        collectMutations(mutated, chain + ',' + currentNode)
      }
    }
  }

  collectMutations(qp, '')

  return { quivers: alreadySeen, chains, timeout: Date.now() - beginTime > maxRuntime }
}

// https://stackoverflow.com/questions/546655/finding-all-cycles-in-a-directed-graph
function findAllCycles (qp, maxCycleLength = 10) {
  var cycles = []

  function collectCyles (qp, originNode, currentNode, visitedNodes, path) {
    if (path.length > maxCycleLength) {
      return
    }
    if (path.length > 0 && currentNode == originNode) {
      cycles.push(path)
    }
    //even once we've completed a cycle, it's possible to create a longer one with a figure-8 shape

    for (var edgeOut of qp.arrowsWithTail[currentNode]) {
      if (!path.includes(edgeOut)) {
        collectCyles(qp, originNode, qp.edges[edgeOut][1], visitedNodes.concat([currentNode]), path.concat([edgeOut]))
      }
    }
  }

  for (var node of qp.nodes) {
    collectCyles(qp, node, node, [], [])
  }

  // cycles = cycles.map(cycle => cycle.sort())

  cycles = cycles.filter(function (cycle, idx) {
    for (var idx2 = 0; idx2 < idx; idx2++) {
      if (JSON.stringify(cycleOrder(cycles[idx])) === JSON.stringify(cycleOrder(cycles[idx2]))) {
        return false
      }
    }
    return true
  })

  return cycles.map(cycle => cycleOrder(cycle))
}

/*
Takes output from findAllCycles, and expands the cycles to have any available combination of self-loops
*/
function extendCyclesWithSelfLoops (cycles, qp, maxCycleLength) {
  var cyclesOut = deepCopy(cycles)

  for (var i = 0; i < cyclesOut.length; i++) {
    var cycle = cyclesOut[i]

    if (cycle.length >= maxCycleLength) {
      continue
    }

    // uncomment to make this non-recursive (limit of 1 self-loop added per cycle)
    // for (var i = 0; i < cycles.length; i++) {
    // var cycle = cycles[i]

    // find if there's a point a self-loop can be inserted here
    for (var e = 0; e < cycle.length; e++) {
      var node = qp.edges[cycle[e]][1]

      for (var e2 = 0; e2 < qp.edges.length; e2++) {
        if (qp.edges[e2][0] === node && qp.edges[e2][1] === node) {
          var cycle2 = deepCopy(cycle)
          cycle2.splice(e + 1, 0, e2)

          // is the new cycle valid?
          // valid if:
          // 1. It doesn't repeat self loops: [6, 8, 8, 7]
          // 2. It doesn't repeat self loops even out of order: [17, 15, 14, 15, 14, 16]
          // (these cycles might still be valid actually, but including them would mean there are infinite possible cycles)
          // 3. It doesn't already exist

          function isSelfLoop (edge) {
            return edge[0] === edge[1]
          }

          var valid = true
          for (var cyclePoint = 1; cyclePoint < cycle2.length; cyclePoint++) {
            if (cycle2[cyclePoint] === cycle2[cyclePoint - 1]) {
              valid = false
              break
            }
            if (isSelfLoop(qp.edges[cycle2[cyclePoint]])) {
              var cp2 = cyclePoint + 1
              while (cycle2[cp2]) {
                if (!isSelfLoop(qp.edges[cycle2[cp2]])) {
                  break
                }
                if (cycle2[cp2] === cycle2[cyclePoint]) {
                  valid = false
                  break
                }
                cp2++
              }
            }
          }

          for (var cout = 0; cout < cyclesOut.length; cout++) {
            if (JSON.stringify(cycleOrder(cyclesOut[cout])) === JSON.stringify(cycleOrder(cycle2))) {
              valid = false
            }
          }

          if (valid) {
            cyclesOut.push(cycle2)
          }
        }
      }
    }
  }

  return cyclesOut.map(cycle => cycleOrder(cycle))
}

function areGraphIsomorphic(q1, q2) {
  if ((q1.nodes.length == q2.nodes.length) && (q1.edges.length == q2.edges.length)) {
    var q1_valences = q1.nodes.map(n => (count(q1.edges.map(e => e[0]), n), count(q1.edges.map(e => e[1]), n)));
    var q2_valences = q2.nodes.map(n => (count(q2.edges.map(e => e[0]), n), count(q2.edges.map(e => e[1]), n)));

    if (q1_valences.sort().join(',') == q2_valences.sort().join(',')) {
      nodeidx = range(0, q1.nodes.length);
      edgeidx = range(0, q1.edges.length);

      m1 = graphAsMat(q1).join(',');
      m2 = graphAsMat(q2);

      var node_permutations = permutations(nodeidx);
      for (var i = 0; i < node_permutations.length; i++) {
	var node_permutation = node_permutations[i];

        var m2_rows = node_permutation.map(i => m2[i]);

        var edge_permutations = permutations(edgeidx);
        for (var j = 0; j < edge_permutations.length; j++) {
	  var edge_permutation = edge_permutations[j];

          var m2_rowscols = m2_rows.map(function(r) {return edge_permutation.map(i => r[i]);});
          if (m1 == m2_rowscols.join(',')) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

function count(l, v) { // return the number of instances of value v in list l
    return l.filter(x => x == v).length;
}

function graphAsMat(quiver) {
    mat_to_return = [];
    for (var n in quiver.nodes) {
        thisrow = [];
        for (var ei = 0; ei < quiver.edges.length; ei++) {
          var e = quiver.edges[ei];
          if (e[0] == n) {thisrow.push(-1);}
          else if (e[1] == n) {thisrow.push(1);}
          else { thisrow.push(0);}
        }
        mat_to_return.push(thisrow);
    }
    return mat_to_return;
}

// https://www.30secondsofcode.org/js/s/permutations
const permutations = arr => {
  if (arr.length <= 2) return arr.length === 2 ? [arr, [arr[1], arr[0]]] : arr;
  return arr.reduce(
    (acc, item, i) =>
      acc.concat(
        permutations([...arr.slice(0, i), ...arr.slice(i + 1)]).map(val => [
          item,
          ...val,
        ])
      ),
    []
  );
};

function range(start, stop, step=1) {
    return Array.from({length: (stop-start)/step}, (_, i) => start+i*step);
}

/*quiver1 = {nodes: [0,1,2,3,4], edges: [[0,1], [1,2], [2,3], [3,4], [4,0], [1,2],]};
quiver2 = {nodes: [0,1,2,3,4], edges: [[0,1], [1,2], [2,3], [3,4], [4,0], [2,3],]};
quiver3 = {nodes: [3,4,0,1,2], edges: [[0,1], [1,2], [2,3], [3,4], [4,0], [2,3],]};
quiver4 = {nodes: [3,4,0,1,2], edges: [[0,1], [1,2], [2,3], [3,4], [4,0], [3,2],]};

areGraphIsomorphic(quiver1, quiver1); // true
areGraphIsomorphic(quiver1, quiver2); // true
areGraphIsomorphic(quiver1, quiver3); // true
areGraphIsomorphic(quiver1, quiver4); // false
*/

function quiverSetsMaybeIsomorphicSimple (setA, setB) {
  // determine if two sets of quivers might be isomorphic to each other

  if (setA.length !== setB.length) {
    return false
  }

  // simple check - edge counts

  var setAStringified = setA.map(q => q.edges.length.toString())
  var setBStringified = setB.map(q => q.edges.length.toString())

  setAStringified.sort()
  setBStringified.sort()

  return setAStringified.join('') === setBStringified.join('')
}

function quiverSetsMaybeIsomorphic (setA, setB) {
  // determine if two sets of quivers might be isomorphic to each other

  if (setA.length !== setB.length) {
    return false
  }

  // simple check - edge counts

  /* var setAStringified = setA.map(q => q.edges.length.toString())
  var setBStringified = setB.map(q => q.edges.length.toString())

  setAStringified.sort()
  setBStringified.sort()

  return setAStringified.join('') === setBStringified.join('') */

  // complex check - edges in and out of each vertex

  function stringifyQuiver (quiver) {
    var edgeCounts = new Array(quiver.nodes.length);
    //can't use array.fill here because it doesn't create unique array instances
    for (var i = 0; i < edgeCounts.length; i++) {
      edgeCounts[i] = [0, 0]
    }

    quiver.edges.forEach(function (e) {
      if (e.hasOwnProperty('from')) {
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

  var setAStringified = setA.map(q => stringifyQuiver(q))
  var setBStringified = setB.map(q => stringifyQuiver(q))

  for (var s of setBStringified) {
    if (setAStringified.includes(s)) {
      setAStringified.splice(setAStringified.indexOf(s), 1)
    } else {
      return false
    }
  }

  return true
}

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

function quiverSetsIsomorphic (setA, setB) {

  if (!quiverSetsMaybeIsomorphicSimple(setA, setB)) {
    return false;
  }
  
  console.log(setA, setB)

  //convert both sets to the right format

  var setAClone = deepCopy(setA)
  var setBClone = deepCopy(setB)

  for (var i = 0; i < setAClone.length; i++) {
    convertQuiver(setAClone[i])
  }
  for (var i = 0; i < setBClone.length; i++) {
    convertQuiver(setBClone[i])
  }

  console.log('after fix', setAClone, setBClone)

  //now, every quiver in set A needs to be isomorphic to some quiver in set B

  for (var i = 0; i < setAClone.length; i++) {
    var found = false;
    for (var i2 = 0; i2 < setBClone.length; i2++) {
      console.log('check', setAClone[i], setBClone[i2])
      if (areGraphIsomorphic(setAClone[i], setBClone[i2])) {
        setBClone.splice(i2, 1)
        found = true
        break;
      }
    }
    if (!found) {
      return false
    }
  }

  return true
}

function isLikelyTerm(qp, cycle) {
  //return arrayEquals(cycle, [1,4,3])// || arrayEquals(cycle, [8,9,10]) || arrayEquals(cycle, [0,2,5]) || arrayEquals(cycle, [2,6,7,3]) || arrayEquals(cycle, [6,8,7]) || arrayEquals(cycle, [0,1,17,16])
 if (cycle.length === 3 && cycle.some(t => qp.edges[t][0] === qp.edges[t][1])) {
    var loops = cycle.filter(t => qp.edges[t][0] === qp.edges[t][1])
    return loops.length === 1
 }
 return false
}

function isLikelyTerm2(qp, cycle) {
  return cycle.length === 3 && !cycle.find(t => qp.edges[t][0] === qp.edges[t][1])
}

function isLikelyTerm3(qp, cycle) {
  if (cycle.length === 4) {
    var tails = {}
    var heads = {}
    cycle.forEach(function(t) {
      const head = qp.edges[t][1]
      const tail = qp.edges[t][0]
      if (tails[tail]) {
        tails[tail]++
      } else {
        tails[tail] = 1
      }
      if (heads[head]) {
        heads[head]++
      } else {
        heads[head] = 1
      }
    })
    for (var node in tails) {
      if (tails[node] === 2 && heads[node] === 2) {
        return true
      }
    }
  }
  return false
}

function cycleIncludesNode(qp, cycle, node) {
  return cycle.some(t => qp.edges[t].includes(node))
}

function potentialRandomSearch (qp, expectedExchangeNum, expectedQuivers = [], maxCycleLength = 5, minPotentialTerms=1, maxPotentialTerms = 100, requiredTerms = [], numberToTest = 5000) {
  // var cyclesWithoutQuadratics = extendCyclesWithSelfLoops(findAllCycles(qp, maxCycleLength), qp).filter(cycle => cycle.length > 2 && cycle.length <= maxCycleLength)
  var cyclesWithoutQuadratics = extendCyclesWithSelfLoops(findAllCycles(qp, maxCycleLength), qp, maxCycleLength).filter(cycle => cycle.length > 2 && cycle.length <= maxCycleLength)

  cyclesWithoutQuadratics = cyclesWithoutQuadratics.filter(c => isLikelyTerm(qp, c) || isLikelyTerm2(qp, c) || isLikelyTerm3(qp, c))

  //specific to 1,2,8
  cyclesWithoutQuadratics = cyclesWithoutQuadratics.filter(cycle => [12, 13, 14, 8].some(node => cycleIncludesNode(qp, cycle, node)))

  /*var termLikeliness = cyclesWithoutQuadratics.map(t => isLikelyTerm(qp, t) || isLikelyTerm2(qp, t))
  var likelyIndexes = termLikeliness.map((i, idx) => idx).filter(idx => termLikeliness[idx] === true)
  var likelyCount = termLikeliness.filter(t => t === true).length*/
  var requiredIndexes = cyclesWithoutQuadratics.map((cycle, idx) => requiredTerms.some(req => arrayEquals(cycleOrder(cycle), cycleOrder(req))) ? idx : undefined).filter(i => i !== undefined)

  console.log(cyclesWithoutQuadratics, requiredIndexes)

  var testedPotentials = []

  var potentialTemplate = cyclesWithoutQuadratics.map(cycle => {
    return [0, cycle.join(',')]
  })
  console.log('testing with terms', cyclesWithoutQuadratics)

  const weightsToTest = [1, 2.5]

  const skipped = 0
  let tested = 0
  let errored = 0
  const erroredResults = []
  // exchange num: count of succeeded
  const exchangeNumBuckets = {}
  const minimalResultsByExchangeNum = {}
  const chainsByExchangeNum = {}
  const maybeMatchingPotentials = []
  const matchingPotentials = []
  const comp = []
  const ruleMatchPotentials = []
  const sizeBuckets = {}

  // limits the terms in the generated potentials to approximately this size

  for (var i = 0; i < numberToTest; i++) {
    var template = deepCopy(potentialTemplate)

    var presetTerms = 0

    /*while (Math.random() < 0.5) {
      var randLikely = likelyIndexes[Math.floor(Math.random() * likelyIndexes.length)]
      template[randLikely][0] = weightsToTest[Math.floor(Math.random() * weightsToTest.length)]
      likeliesSet++
    }*/
    /*for (var li = 0; li < likelyIndexes.length; li++) {
      if (Math.random() < Math.pow(0.5, likeliesSet + 1)) {
        template[likelyIndexes[li]][0] = weightsToTest[Math.floor(Math.random() * weightsToTest.length)]
        likeliesSet++
      }
    }*/
    for (var li = 0; li < requiredIndexes.length; li++) {
      //if (Math.random() < 0.9) {
       // template[requiredIndexes[li]][0] = weightsToTest[Math.floor(Math.random() * weightsToTest.length)]
       template[requiredIndexes[li]][0] = 1
        presetTerms++
     // }
    }

  /*  for (var li = 0; li < likelyIndexes.length; li++) {
      if (Math.random() < (2 / likelyCount)) {
        template[likelyIndexes[li]][0] = weightsToTest[Math.floor(Math.random() * weightsToTest.length)]
        presetTerms++
      }
    }*/

    const thisPotentialSize = presetTerms + Math.round(Math.random() * (maxPotentialTerms - presetTerms))

    const thisPotentialFactor = Math.min(1, (thisPotentialSize - presetTerms) / cyclesWithoutQuadratics.length)

    // this makes the generated potentials linearly distributed with respect to their size
    /*var thisPotentialFactor = Math.random() * potentialAdjustFactor

    const likelyFactor = Math.min(thisPotentialFactor * 10, 1)
    const unlikelyFactor = ((thisPotentialFactor * cyclesWithoutQuadratics.length) - (likelyFactor * likelyCount)) / (cyclesWithoutQuadratics.length - likelyCount)

    for (var t = 0; t < template.length; t++) {
      if ((termLikeliness[t] && Math.random() < thisPotentialFactor * likelyFactor) || 
         (!termLikeliness[t] && Math.random() < unlikelyFactor)) {
        template[t][0] = weightsToTest[Math.floor(Math.random() * weightsToTest.length)]
      }
    }*/

    for (var t = 0; t < template.length; t++) {
      if (Math.random() < thisPotentialFactor) {
        template[t][0] = weightsToTest[Math.floor(Math.random() * weightsToTest.length)]
      }
    }

    /*if (thisPotentialSize === 6 && likeliesSet === 6) {
      console.log(template)
    }*/

    // uncomment to skip duplicate potentials
    // var templateStr = JSON.stringify(template)
    // if(testedPotentials.includes(templateStr)) {
    //     continue;
    // }
    // testedPotentials.push(templateStr);

    var qpt = deepCopy(qp)
    var constructedPotential = template.filter(t => t[0] !== 0)

    if (constructedPotential.length < minPotentialTerms || constructedPotential.length > maxPotentialTerms) {
      i--
      continue
    }

    if (sizeBuckets[constructedPotential.length]) {
      sizeBuckets[constructedPotential.length]++
    } else {
      sizeBuckets[constructedPotential.length] = 1
    }

/*     if (sizeBuckets[constructedPotential.length]) {
      sizeBuckets[constructedPotential.length]++
    } else {
      sizeBuckets[constructedPotential.length] = 1
    }
 */
    qpt.potential = constructedPotential
    try {
      var exchangeNumResult = getAllMutationsForQP(qpt, expectedExchangeNum + 1)
      var exchangeNum = exchangeNumResult.quivers.length

      if (constructedPotential.some(t => t[1] === '0,2,5') && constructedPotential.some(t => t[1] === '1,4,3') && constructedPotential.some(t => t[1] === '8,9,10') && constructedPotential.some(t => t[1] === '2,6,7,3') && constructedPotential.some(t => t[1] === '0,1,17,16') && constructedPotential.some(t => t[1] === '6,8,7')) {
        //ruleMatchPotentials.push(constructedPotential)
      }

      if (exchangeNum === expectedExchangeNum) {
        if (quiverSetsMaybeIsomorphic(exchangeNumResult.quivers.map(qp => JSON.parse(qp)), expectedQuivers)) {
          console.log('likely', JSON.stringify(constructedPotential))
          maybeMatchingPotentials.push(constructedPotential)
        }
        /*if (quiverSetsIsomorphic(exchangeNumResult.quivers, expectedQuivers)) {
          matchingPotentials.push(constructedPotential)
        }*/
      }

      // if (chainsByExchangeNum[exchangeNum]) {
      //    chainsByExchangeNum[exchangeNum].push(exchangeNumResult.chains)
      // } else {
      //    chainsByExchangeNum[exchangeNum] = [exchangeNumResult.chains]
      // }

      if (exchangeNumBuckets[exchangeNum]) {
        exchangeNumBuckets[exchangeNum]++
      } else {
        exchangeNumBuckets[exchangeNum] = 1
      }

      if (!minimalResultsByExchangeNum[exchangeNum]) {
        minimalResultsByExchangeNum[exchangeNum] = [constructedPotential]
      } else if (constructedPotential.length < minimalResultsByExchangeNum[exchangeNum][0].length) {
        minimalResultsByExchangeNum[exchangeNum] = [constructedPotential]
      } else if (constructedPotential.length === minimalResultsByExchangeNum[exchangeNum][0].length) {
        minimalResultsByExchangeNum[exchangeNum].push(constructedPotential)
      }
    } catch (e) {
      //   console.log(e)
      errored++
      //  erroredResults.push(constructedPotential)
    }

    tested++
    // console.log(tested);
    if (tested % 1000 === 0) {
      console.log(tested / numberToTest, '%', errored)
    }
  }

  for (var key in minimalResultsByExchangeNum) {
    // deduplicate the minimal results
    var stringifiedResults = minimalResultsByExchangeNum[key].map(i => JSON.stringify(i))
    minimalResultsByExchangeNum[key] = minimalResultsByExchangeNum[key].filter((item, idx) => stringifiedResults.indexOf(JSON.stringify(item)) === idx)
  }

  return {
    stats: {
      tested,
      skipped,
      errored
    },
    // erroredResults,
    exchangeNumBuckets,
    minimalResultsByExchangeNum,
    chainsByExchangeNum,
    maybeMatchingPotentials,
    matchingPotentials,
    comp,
    ruleMatchPotentials,
    sizeBuckets
  }
}

function diffQPEdges(qpA, qpB) {
  var notInA = deepCopy(qpB.edges)
  var notInB = deepCopy(qpA.edges)

  for (var i = 0; i < qpA.edges.length; i++) {
    var idx = notInA.findIndex(edg => arrayEquals(edg, qpA.edges[i]))
    if (idx !== -1) {
      notInA.splice(idx, 1)
    }
  }

  for (var i = 0; i < qpB.edges.length; i++) {
    var idx = notInB.findIndex(edg => arrayEquals(edg, qpB.edges[i]))
    if (idx !== -1) {
      notInB.splice(idx, 1)
    }
  }

  return {notInA, notInB}
}

function eliminateType1(oldQP, extraEdges, lastMut) {
  /*
  case 1: a pair of shortcut edges between two nodes adjacent to the mutated one */

  var adjacentNodes = oldQP.edges.filter(edg => edg[0] === lastMut || edg[1] === lastMut).map(edg => edg[0] === lastMut ? edg[1] : edg[0])
  adjacentNodes = adjacentNodes.filter((i, idx) => adjacentNodes.indexOf(i) === idx)
  

  var adjacentPairs = permutations(adjacentNodes).map(p => p.slice(0, 2))
  //this is inefficient (includes duplicate pairs) but works for now as long as node counts are small

  //remove adjacent node pairs that connect to each other - that is type 2

  adjacentPairs = adjacentPairs.filter(function(pair) {
    return !oldQP.edges.some(edg => arrayEquals(edg, pair))
  })

  for(var i = 0; i < adjacentPairs.length; i++) {
    var thisPair = adjacentPairs[i]
    var idx1 = extraEdges.findIndex(edg => arrayEquals(edg, thisPair))
    var idx2 =  extraEdges.findIndex(edg => arrayEquals(edg, deepCopy(thisPair).reverse()))
    if (idx1 !== -1 && idx2 !== -1) {
      extraEdges.splice(idx1, 1)
      //need to recalc because the first splice changed the index
      idx2 =  extraEdges.findIndex(edg => arrayEquals(edg, deepCopy(thisPair).reverse()))
      extraEdges.splice(idx2, 1)
      return [
          [
          1,
          [
            [thisPair[0], lastMut],
            [lastMut, thisPair[1]],
            [thisPair[1], lastMut],
            [lastMut, thisPair[0]]
          ].map(function(edg) {
              return oldQP.edges.findIndex(t => arrayEquals(t, edg))
            }).join(",")
        ]
    ]
    }
  }
}

function eliminateType2(oldQP, extraEdges, lastMut) {
  /*
  case 2: two adjacent nodes that are connected and shouldn't be
   */

  var adjacentNodes = oldQP.edges.filter(edg => edg[0] === lastMut || edg[1] === lastMut).map(edg => edg[0] === lastMut ? edg[1] : edg[0])
  adjacentNodes = adjacentNodes.filter((i, idx) => adjacentNodes.indexOf(i) === idx)
  
  var adjacentPairs = permutations(adjacentNodes).map(p => p.slice(0, 2))
  //this is inefficient (includes duplicate pairs) but works for now as long as node counts are small

  //remove duplicates
  adjacentPairs = adjacentPairs.filter((pair, i) => adjacentPairs.findIndex(p => arrayEquals(p, pair)) === i)

  //the pair connects to each other with extra edges

  for (var p = 0; p < adjacentPairs.length; p++) {
    var pair = adjacentPairs[p]
    if (extraEdges.filter(edg => arrayEquals(edg, pair)).length === 2
    && extraEdges.filter(edg => arrayEquals(edg, deepCopy(pair).reverse())).length === 2) {
      //we found something
      //remove the edges
      for (var i = 0; i < extraEdges.length; i++) {
        if (arrayEquals(extraEdges[i], pair) || arrayEquals(extraEdges[i], deepCopy(pair).reverse())) {
          extraEdges.splice(i, 1);
          i--;
        }
      }
      return [
        [
          1,
          [
            [lastMut, pair[0]],
            [pair[0], pair[1]],
            [pair[1], lastMut]
          ].map(function(edg) {
            return oldQP.edges.findIndex(t => arrayEquals(t, edg))
          }).join(",")
        ],
        [
          1,
          [
            [lastMut, pair[1]],
            [pair[1], pair[0]],
            [pair[0], lastMut]
          ].map(function(edg) {
            return oldQP.edges.findIndex(t => arrayEquals(t, edg))
          }).join(",")
        ]
    ]
    }
  }
}

function eliminateType3(oldQP, extraEdges, lastMut) {
  //double self loop adjacent
  var adjacentNodes = oldQP.edges.filter(edg => edg[0] === lastMut || edg[1] === lastMut).map(edg => edg[0] === lastMut ? edg[1] : edg[0])
  adjacentNodes = adjacentNodes.filter((i, idx) => adjacentNodes.indexOf(i) === idx)

  var adjacentSelfLoops = extraEdges.filter(edg => edg[0] === edg[1] && adjacentNodes.includes(edg[0]))

  for(var i = 0; i < adjacentNodes.length; i++) {
    if (adjacentSelfLoops.filter(edg => edg[0] === adjacentNodes[i]).length === 2) {
      //remove two self loops
      extraEdges.splice(extraEdges.findIndex(edg => edg[0] === adjacentNodes[i] && edg[1] === adjacentNodes[i]), 1)
      extraEdges.splice(extraEdges.findIndex(edg => edg[0] === adjacentNodes[i] && edg[1] === adjacentNodes[i]), 1)

      return [
        [
          1, 
          [
            [lastMut, adjacentNodes[i]],
            [adjacentNodes[i], adjacentNodes[i]],
            [adjacentNodes[i], lastMut]
          ].map(function(edg) {
            return oldQP.edges.findIndex(t => arrayEquals(t, edg))
          }).join(",")
        ]
      ]
    }
  }
}

function potentialStructuredSearch(expectedQuivers, mutationSequences) {
  if (mutationSequences[0].length !== 0) {
    throw new Error("first QP must be the base")
  }

  expectedQuivers = remapQPNodes(expectedQuivers).map(q => convertQuiver(q)).map(base => makeQP(base.edges, base.nodes, base.frozenNodes, base.potential, 'fromThing'))

  var expectedWithSequences = expectedQuivers.slice(1).map((q, i) => {return {q, seq: mutationSequences[i + 1]}})

  expectedWithSequences = expectedWithSequences.sort(function(a, b) {
    return a.seq.length - b.seq.length
  })

  //console.log(expectedWithSequences)

  //build a tree of the mutation sequences

  /*var mutationTree = {}
  mutationSequences.forEach(function(seq) {
    var ref = mutationTree
    seq.forEach(function(node) {
      if (!ref[node]) {
        ref[node] = {}
      }
      ref = ref[node]
    })
    ref.done = true
  })*/

  var existingTerms = []

  expectedWithSequences.forEach(function(obj) {
    var qp = deepCopy(expectedQuivers[0])
    var lastQP = qp
    qp.potential = deepCopy(existingTerms)
   // console.log(qp.potential)

    obj.seq.forEach(function(node) {
     // console.log('mutate', obj.seq, node)
      lastQP = qp
      qp = mutateQP(node, deepCopy(qp))
    })

    //now find the extra edges

    var diff = diffQPEdges(qp, obj.q)

    //now add terms

    var extraEdges = diff.notInB
   // console.log('begin with', deepCopy(existingTerms))

    var newTerms = []

    while (extraEdges.length > 0) {
    //  console.log('at iter', deepCopy(extraEdges))
      var t1 = eliminateType1(lastQP, extraEdges, obj.seq[obj.seq.length - 1])

      if (t1) {
    //    console.log('type 1 elimination', t1)
        newTerms = newTerms.concat(t1)
        continue
      }

      var t2 = eliminateType2(lastQP, extraEdges, obj.seq[obj.seq.length - 1])

      if (t2) {
     //   console.log('type 2 elimination', t2)
        newTerms = newTerms.concat(t2)
        continue
      }

      var t3 = eliminateType3(lastQP, extraEdges, obj.seq[obj.seq.length - 1])

      if (t3) {
     //   console.log('type 3 elimination', t3)
        newTerms = newTerms.concat(t3)
        continue
      }

      //we weren't able to eliminate anything
      throw new Error("elimination failed")
    }

   // console.log("ended with", deepCopy(existingTerms), deepCopy(newTerms))

    //now mutate the lastQP with the terms back to base

    var newQP = deepCopy(lastQP)
    newQP.potential = deepCopy(lastQP.potential).concat(newTerms)

    obj.seq.slice(0, -1).reverse().forEach(function(node) {
    //  console.log('mutate back', node)
    //  console.log(tryBuildReplacementPotentials(deepCopy(newQP)))
      newQP = mutateQP(node, newQP)
    })
    //console.log(tryBuildReplacementPotentials(deepCopy(newQP)))

    //console.log('done', obj.seq, newQP)

    //now mutate back to base

    //now we have a new potential
    //remap the edge IDs back to original

    var remappedPotential = newQP.potential.map(function(term) {
      var newTerm = term[1].split(",").map(function(edg) {
        var oldEdg = newQP.edges[edg];
        var newEdg = expectedQuivers[0].edges.findIndex(n => arrayEquals(n, oldEdg))
        return newEdg
      }).join(",")
      return [
        term[0],
        newTerm
      ]
    })

    //console.log(deepCopy(newQP.potential), deepCopy(remappedPotential))

    existingTerms = deepCopy(remappedPotential)
  })

  //find alternatives?

  var altTestQP = deepCopy(expectedQuivers[0])
  altTestQP.potential = deepCopy(existingTerms)
  console.log(tryBuildReplacementPotentials(altTestQP))

  return existingTerms

 // console.log('final result', existingTerms)
}

/*
Used by potential structured search
Given a set of (randomly-generated) mutation chains and the resulting potential, simplify it by elimating the terms that don't do anything
(because there is a loop at that node, so that mutation step has no effect)
*/
function simplifyMutationChains(baseQP, potential, mutationChains) {
    baseQP = convertQuiver(deepCopy(baseQP))
    baseQP = makeQP(baseQP.edges, baseQP.nodes, baseQP.frozenNodes, baseQP.potential, 'fromThing')
    baseQP.potential = potential
    return mutationChains.map(function(chain) {
      var tqp = deepCopy(baseQP)
      var outChain = []
      chain.forEach(function(chainElement) {
        if (tqp.canMutate[chainElement]) {
          outChain.push(chainElement)
          tqp = mutateQP(chainElement, tqp)
        }
      })
      return outChain
    })
}

https://stackoverflow.com/questions/6274339/how-can-i-shuffle-an-array
function shuffleArr(a) {
  var j, x, i;
  for (i = a.length - 1; i > 0; i--) {
      j = Math.floor(Math.random() * (i + 1));
      x = a[i];
      a[i] = a[j];
      a[j] = x;
  }
  return a;
}

function findPossibleMutationNodes(qp1, qp2) {
  var qps = remapQPNodes(deepCopy([qp1, qp2]))
  .map(q => convertQuiver(q)).map(base => makeQP(base.edges, base.nodes, base.frozenNodes, base.potential, 'fromThing'))

  var edgeDiff = diffQPEdges(qps[0], qps[1])

  var possibleNodes = null
  edgeDiff.notInB.forEach(function(edg) {
    var adjacent = qps[0].nodes.filter(function(node) {
      return node === edg[0] || node === edg[1] || qps[0].edges.some(function(otherEdge) {
        return (otherEdge[0] === node && edg.includes(otherEdge[1]))
        || (otherEdge[1] === node && edg.includes(otherEdge[0]))
      })
    })
    if (possibleNodes === null) {
      possibleNodes = adjacent
    } else {
      possibleNodes = possibleNodes.filter(n => adjacent.includes(n))
    }
  })

  edgeDiff.notInA.forEach(function(edg) {
    var adjacent = qps[1].nodes.filter(function(node) {
      return node === edg[0] || node === edg[1] || qps[1].edges.some(function(otherEdge) {
        return (otherEdge[0] === node && edg.includes(otherEdge[1]))
        || (otherEdge[1] === node && edg.includes(otherEdge[0]))
      })
    })
    if (possibleNodes === null) {
      possibleNodes = adjacent
    } else {
      possibleNodes = possibleNodes.filter(n => adjacent.includes(n))
    }
  })
  if (possibleNodes === null) {
    console.warn('was passed two identical quivers to compare')
  }
  return possibleNodes || []
}

function findMutationChainsForQPSet(allQPs) {
  var baseMutationChains = allQPs.map(function(qp, i) {
    if (i === 0) {
      return []
    } else {
      return null;
    }
  })

    outer: while (true) {      
       for(var i = 0; i < baseMutationChains.length; i++) {
        if (baseMutationChains[i] === null) {
          for (var i2 = 0; i2 < baseMutationChains.length; i2++) {
            if (i !== i2 && baseMutationChains[i2] !== null) {
              var possible = findPossibleMutationNodes(allQPs[i2], allQPs[i])
              if (possible.length === 1) {
                baseMutationChains[i] = deepCopy(baseMutationChains[i2]).concat(possible)
                continue outer;
              }
            }
          }
        }
      }
      //wasn't able to find anything
      break outer;
    }
    return baseMutationChains
}


function potentialStructuredRandomSearch(allQPs, iter=10000, maxDepth = 3) {
  //we know that any node mutable in the base QP must generate another QP in the result
  var requiredMutationNodes = []
  allQPs[0].nodes.forEach(function(node) {
    var hasLoop = allQPs[0].edges.some(edg => edg.from === node.id && edg.to === node.id)
    if (!hasLoop) {
      requiredMutationNodes.push(parseInt(node.id))
    }
  })

  var baseMutationChains = findMutationChainsForQPSet(allQPs)

  //can't generate anything random
  if (!baseMutationChains.some(i => i === null)) {
    iter = 1;
  }

  //TODO we should be able to match all of the required chains to a QP just based on which edges change
  //keeping in mind that edges need to be remapped
  //for a single edge to be added or removed, an adjacent node or the third in a triangle must have been mutated
  //take the intersection of all of these
  //can we repeat this for the remaining quivers?
  //and if we do base -> a 2+ quiver with this method, do we detect that we can't reach it?
  outer: for (var i = 0; i < iter; i++) {
    var mutationChains = deepCopy(baseMutationChains)
     if (i % 1000 === 0) {
      console.log(i)
    }

    for (var c = 0; c < mutationChains.length; c++) {
      if (mutationChains[c] === null) {
        var arr = new Array(Math.round(Math.random() * (maxDepth - 1))).fill(0)
        arr = arr.map(i => Math.floor(Math.random() * allQPs[0].nodes.length))
        //any mutation chain must start with a node that's initially mutable
        arr.unshift(requiredMutationNodes[Math.floor(Math.random() * requiredMutationNodes.length)])
        mutationChains[c] = arr
      }
    }
      
    try {
      var result = potentialStructuredSearch(allQPs, mutationChains);
      var simplifiedChains = simplifyMutationChains(allQPs[0], result, mutationChains)
      return [result, simplifiedChains, i]
      break outer;
    } catch (e) {}
  }
}

function potentialStructuredTest(max=100) {
  var results = {
    failedTriangulation: [],
    exchangeNumTooBig: [],
    failedGenerate: [],
    failedCheck: [],
    successes: [],
  }
  for (var a = 1; a <= max; a++) {
    for(var b= 1; b <= max; b++) {
      for (var c = 1; c <= max; c++) {
        console.log(a, b, c)
        var succeeded = false;
        var trials = 0;
        /*
        It's possible for potentialStructuredRandomSearch to produce a potential and sequence of mutation chains that isn't right
        I believe the issue is that the mutation chains are generated blindly, and sometimes an impossible chain is generated (because there's a loop at the location that should be mutated) - also the final potentail affects which chains are possible
        So we perform an extra step to verify that the result is correct, and try again if not
        */
        while (!succeeded && trials < 50) {
          trials++
          try {
          var r = a + b + c
          var tri = triangulation(r, a, b, c);
          var data = allUniqueTriangulations(tri, getBoundaryEdges(r, tri[0], tri[1]))
          const cs = data.coordinates.map(y => y.map(z=>parseFloat(z)));
          data = data.triangulations.map(x => JSON.parse(x));
          data = data.map(t => t.map(e => JSON.parse(e)));
          data = data.map(x => QPFromTriangulation([x, cs]));
          //TODO remove
          if (data.length > 500) {
            console.warn('skipping ', [r, a, b, c].join(",") + " because the exchange number is too big")
            results.exchangeNumTooBig.push([r, a, b, c])
            break
          }
          } catch (e) {
            results.failedTriangulation.push([r, a, b, c])
            break
          }

          var result = potentialStructuredRandomSearch(data)
          if (!result) {
            results.failedGenerate.push([r, a, b, c])
            break;
          } 
          var verification = doPartialComparison(result[0], data)
          
           if (verification[2] === true) {
            results.successes.push({
              input: [r, a, b, c],
              data: data,
              output: result,
              verification
            })
            succeeded = true
          }
        }
        if (!succeeded) {
          results.failedCheck.push([r, a, b, c])
        }
      }
    }
  }
  return results;                        
}

function potentialIncludesEveryVertex (qp, potential) {
  const nodesInPotential = []
  potential.forEach(function (term) {
    if (term[0] !== 0) {
      term[1].split(',').map(i => parseInt(i)).forEach(function (edge) {
        qp.edges[edge].forEach(function (node) {
          if (!nodesInPotential.includes(node)) {
            nodesInPotential.push(node)
          }
        })
      })
    }
  })
  return nodesInPotential.length === qp.nodes.length
}

function isPotentialSubsetOf (potA, potB) {
  for (var a = 0; a < potA.length; a++) {
    var existsInB = false
    for (var b = 0; b < potB.length; b++) {
      if (potB[b][1] === potA[a][1] && (potA[a][0] === 0 || potB[b][0] !== 0)) {
        existsInB = true
      }
    }
    if (!existsInB) {
      return false
    }
  }
  return true
}

function pathDerivative (thisPotential, edgeIndex, fmt = 'string') {
  if (fmt == 'string') {
    var tp = thisPotential.filter(function (tc) {
      return tc[1].split(',').indexOf(edgeIndex.toString()) >= 0
    }).map(function (termcoef) {
      const allTerms = termcoef[1].split(',').map(y => parseInt(y))
      const ati = allTerms.indexOf(edgeIndex)
      return [termcoef[0], allTerms.slice(ati + 1).concat(...allTerms.slice(0, ati))]
    })
  } else {
    var tp = thisPotential.filter(function (tc) {
      return tc[1].includes(edgeIndex)
    }).map(function (termcoef) {
      const allTerms = termcoef[1]
      const ati = allTerms.indexOf(edgeIndex)
      return [termcoef[0], allTerms.slice(ati + 1).concat(...allTerms.slice(0, ati))]
    })
  }
  if (tp != null) { return tp } else { return [] }
}

function randInt (range) {
  const sign = Math.random() > 0.5 ? 1 : -1
  return sign * Math.floor(Math.random() * range)
}

function randomPotential (ns, es, coefficient_range = 100) {
  theseNodes = ns.getIds().map(x => parseInt(x))
  theseEdges = es.getIds().map(x => [parseInt(es.get(x).from), parseInt(es.get(x).to)])
  const thisQP = makeQP(theseEdges, theseNodes, [0], [[0, '0']], 'fromThing')
  const allTriples = allThreeCycles(thisQP)

  var currentlyMetNodes = []
  var metAllNodes = false
  var cycles = []
  var metCycles = []
  do {
    let randomCycle = randInt(allTriples.length)
    randomCycle = randomCycle >= 0 ? randomCycle : -randomCycle
    const theCycle = allTriples[randomCycle]

    if (!metCycles.includes(randomCycle)) {
      for (let tci = 0; tci < 3; tci++) {
        const e = theseEdges[theCycle[tci]]
        if (!currentlyMetNodes.includes(e[0])) {
          currentlyMetNodes.push(e[0])
        }
        if (!currentlyMetNodes.includes(e[1])) {
          currentlyMetNodes.push(e[1])
        }
      }
      metCycles.push(randomCycle)
      cycles.push(theCycle.toString())
    }

    if ((currentlyMetNodes.length >= ns.length) || (metCycles.length >= allTriples.length)) {
      metAllNodes = true
    }
  } while (!metAllNodes)

  return cycles.map(x => [randInt(coefficient_range), x])
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

function reduce (QP) {
  // remove extraneous commas from potential
    var thePotential = QP.potential.map(function(x) {
      return [Number(x[0]), x[1].split(',').filter(i => i !== null).map(z => parseInt(z))]
    })
    .filter(termcoef => (Math.abs(termcoef[0]) > 0 && termcoef[1].length > 0))

  // extract the terms that show up as a singleton after taking a path derivative. (These terms are equivalent to 0)
  var allPathDerivatives = range(0, QP.edges.length).map(e => pathDerivative(thePotential, e, fmt='list')).filter(pd => pd.length > 0)
  var zeroTerms = allPathDerivatives.filter(pdterm => (pdterm.length == 1)).map(pdt => pdt[0][1])
  zeroTerms = zeroTerms.filter(t => (t.length == 1)).map(t => t[0])

  // now remove all terms that contain an edge equivalent to 0
  if (zeroTerms.length > 0) {
    thePotential = thePotential.filter(function (termcoef) {
      return termcoef[1].some(t => zeroTerms.indexOf(t) < 0)
    })
  }

  function termsContainEdge (terms, edge) {
    return terms.some(x => x[1].includes(edge))
  }
  var squareTerms = thePotential.filter(x => x[1].length == 2).map(y => y[1])
  var squareCoefs = thePotential.filter(x => x[1].length == 2).map(y => y[0])
  var edgesToRemove = unique(squareTerms.flat())

  if (squareTerms.length > 0) {
    // create a 'lookup dictionary' of replacement terms as follows:
    // reduceDict[edge] = [...] where [...] is either:
    //     1. e1 itself (if e1 is not an edge that occurs in a quadratic term, then we won't try to remove it)
    //     2. a term that is equivalent to e1 in the jacobian algebra (if we're removing e1)
    var reduceDict = range(0, QP.edges.length).map(x => [[1, [x]]])
    var alreadySeen = range(0, QP.edges.length).map(x => false)
    for (let ti = 0; ti < squareTerms.length; ti++) {
	    const t = squareTerms[ti]
	    const c = squareCoefs[ti]
	    const e1 = t[0]; const e2 = t[1]

      if (!alreadySeen[e1]) {
	        alreadySeen[e1] = true
        reduceDict[e1] = pathDerivative(thePotential, e2, fmt = 'list').map(
          function (x) {
            if (!((x[1].filter(y => y != null).length < 2) && x[1].includes(e1))) {
              return [-x[0] / c, x[1]]
            }
          }).filter(y => (y != null))

        // double check that this new term term is not of the form A = AX + B
        // (i.e. the replacement terms for edge A does not contain a term with A in it)
        if (termsContainEdge(reduceDict[e1], e1)) {
          edgesToRemove.splice(edgesToRemove.indexOf(e1), 1)
          reduceDict[e1] = [[1, [e1]]]
        }
      }
      if (!alreadySeen[e2]) {
	        alreadySeen[e2] = true
        reduceDict[e2] = pathDerivative(thePotential, e1, fmt = 'list').map(
          function (x) {
            if (!((x[1].filter(y => y != null).length < 2) && x[1].includes(e2))) {
              return [-x[0] / c, x[1]]
            }
          }).filter(y => (y != null))

        if (termsContainEdge(reduceDict[e2], e2)) {
          edgesToRemove.splice(edgesToRemove.indexOf(e2), 1)
          reduceDict[e2] = [[1, [e2]]]
        }
      }
    }

    const maxTermsForE = 250
    const maxEdgesPerTerm = 250

    // now recursively replace terms in reduceDict's values so that every edge-to-be-removed
    // is replaced with terms that don't contain edges-to-be-removed.

    var failedToReplace = []
    for (let edgeToRemovei = 0; edgeToRemovei < edgesToRemove.length; edgeToRemovei++) {
      const e = edgesToRemove[edgeToRemovei]
      var termsForE = reduceDict[e].map(x => [x[0], [...x[1]]])
      var foundReplacement = true
      var ctr = 0

      do {
        if (termsForE.length > maxTermsForE) {
          throw new Error('reduction timeout')
        }
        // stopping criteria
        ctr += 1; foundReplacement = true

        // placeholder for holding non-edgesToRemove lookup values for edge e
        var altTermsForE = []
        for (let cti = 0; cti < termsForE.length; cti++) {
          const currentTerm = termsForE[cti]

          if (currentTerm.length > 0) {
            var altTerm = [[currentTerm[0], currentTerm[1]]]

            // check if any of the terms in e's replacement
            // terms also contains one of the edges to remove
            if (currentTerm[1] != null) {
              if (currentTerm[1].some(x => (edgesToRemove.includes(x)))) {
                foundReplacement = false

                // if so, then we need to replace that term
                var newTerm = [[1, []]]
                for (let ttt = 0; ttt < currentTerm[1].length; ttt++) {
                  const tt = currentTerm[1][ttt]
                  if (edgesToRemove.includes(tt) && !failedToReplace.includes(tt)) {
                    var nt = []
                    for (let rdi = 0; rdi < reduceDict[tt].length; rdi++) {
                      const rd = reduceDict[tt][rdi]
                      for (let nt1i = 0; nt1i < newTerm.length; nt1i++) {
                        if (currentTerm[1].length > maxEdgesPerTerm || newTerm.length > maxEdgesPerTerm) {
                          throw new Error('reduction timeout - 2')
                        }

                        const nt1 = newTerm[nt1i]
                        var nt11 = rd[1]
                        if (nt1[1].length > 0) {
                          nt11 = nt1[1].concat(rd[1])
                        }
                        nt.push([parseFloat(nt1[0]) * parseFloat(rd[0]), nt11])
                      }
                    }
                    if (nt.length > 0) { newTerm = nt }
                  } else {
                    newTerm = newTerm.map(
                      function (x) {
                        if (x[1].length > 0) {
                          return [x[0], x[1].concat(tt)]
                        } else {
                          return [x[0], [tt]]
                        }
                      })
                  }
                }
                if (newTerm[0][1].length > 0) {
                  altTerm = newTerm.map(x => [x[0] * currentTerm[0], x[1]])
                }
              }
            }
            altTermsForE.push(...altTerm)
          }
        }
        if (termsContainEdge(termsForE, e)) {
          foundReplacement = false
          ctr = edgesToRemove.length + 1
          termsForE = [[1, [e]]]
          altTermsForE = [[1, [e]]]
        } else {
          reduceDict[e] = simpleArrayClone(termsForE)
        	    termsForE = simpleArrayClone(altTermsForE)
        }
      } while (!foundReplacement && (ctr < edgesToRemove.length))
	    reduceDict[e] = termsForE

      if (!foundReplacement) {
        failedToReplace.push(e)
        reduceDict[e] = [[1, [e]]]
      }
    }

    for (let e2r = 0; e2r < failedToReplace.length; e2r++) {
      const edgeToRemove = failedToReplace[e2r]
      if (edgesToRemove.indexOf(edgeToRemove) >= 0) {
        edgesToRemove.splice(edgesToRemove.indexOf(edgeToRemove), 1)
      }
    }

    const maxTermLength = 250

    // reduce the potential by replacing each of the terms with its
    // image in the replacement dictionary.
    var wPrime = []
    for (let tci = 0; tci < thePotential.length; tci++) {
      const coef = thePotential[tci][0]
      const term = thePotential[tci][1]

      var newTerm = [[coef, []]]
      for (let edgeInTermi = 0; edgeInTermi < term.length; edgeInTermi++) {
        const edgeInTerm = term[edgeInTermi]

        var thisLevel = []
        for (let nt1i = 0; nt1i < newTerm.length; nt1i++) {
          if (term.length > maxTermLength || newTerm.length > maxTermLength) {
            throw new Error('term too big')
          }
          var nt1 = newTerm[nt1i]

          for (let rdi = 0; rdi < reduceDict[edgeInTerm].length; rdi++) {
            var rd = reduceDict[edgeInTerm][rdi]
            thisLevel.push([parseFloat(nt1[0]) * parseFloat(rd[0]),
              nt1[1].concat(rd[1])])
          }
        }
        if (thisLevel.length > 0) { newTerm = simpleArrayClone(thisLevel) }
      }
      wPrime.push(...newTerm)
    }
    wPrime = wPrime.filter(y => (y[0] != 0)).map(function (x) { return [x[0], x[1].toString()] })
    wPrime = combineLikeTermsInPotential(wPrime)

    return removeEdges(edgesToRemove, QP, altPotential = wPrime)
  } else {
    return deepCopy(QP)
  }
}

function range (start, stop, step = 1) { // trying to be as python-ish as I can
  return Array.from({ length: (stop - start) / step }, (_, i) => start + (i * step))
}

function removeEdges (edgeIndices, QP, altPotential = 'None') {
  var edgesToKeep = [...Array(QP.edges.length).keys()].map(
    function (x) {
      if (edgeIndices.includes(parseInt(x))) {
        return -1
      } else {
        return parseInt(x)
      }
    })
  var edgeIndexLookup = edgesToKeep.filter(x => x >= 0)
  var edgeIndexLookupBackwards = edgesToKeep.map(function (x) { if (x >= 0) { return edgeIndexLookup.indexOf(x) } })

  // re-index the edges to delete those that are no longer included
  var newEdges = edgeIndexLookup.map(x => QP.edges[x])

  // update the terms in the potential to use the new edge indexing
  var newPotential = altPotential
  if (newPotential == 'None') {
    newPotential = QP.potential.map(x => [...x])
  }
  newPotential = newPotential.map(
    function (x) {
	    if (x[1].split(',').filter(y => edgeIndices.includes(parseInt(y))).length < 1) {
        const y = x[1].split(',').map(y => edgeIndexLookupBackwards[parseInt(y)])
	        return [parseFloat(x[0]), y.filter(x => x != null).toString()]
      }
    }).filter(x => x != null)
  return makeQP(newEdges, QP.nodes, QP.frozenNodes, newPotential, inputType = 'fromQP')
}

if (typeof module !== 'undefined') {
  module.exports = {makeQP, stringifyQP, mutateQP, potentialRandomSearch}
}

//find 3-cycles that don't have an inverse, and could therefore be replaced by a 4-cycle (assuming it's at a mutable node, which this doens't check)
function tryBuildReplacementPotentials(qp) {
  var suggestedReplacements = []

  var threeTerms = qp.potential.filter(function(term) {
    var edges = term[1].split(",").map(e => parseInt(e))
    var edges_nonloops = edges.filter(edg => qp.edges[edg][0] !== qp.edges[edg][1])
    return edges.length === 3 && edges_nonloops.length === 3
  })

  threeTerms = threeTerms.filter(function(term) {
    var edges = term[1].split(",").map(e => parseInt(e))
    //try to reverse each edge
    var edgesReversed = edges.map(function(edg) {
      var thisEdge = qp.edges[edg]
      var otherEdge = qp.edges.findIndex(other => other[1] === thisEdge[0] && other[0] === thisEdge[1])
      return otherEdge
    })
    //couldn't reverse all the edges
    if (edgesReversed.includes(undefined)) {
      return true
    }
    //see if the reverse term exists
    var reverseTermExists = threeTerms.some(function(otherTerm) {
      var other = deepCopy(otherTerm[1]).split(",").sort().join(",")
      var thisTerm = edgesReversed.map(t => t.toString()).sort().join(",")
      return other === thisTerm
    })

    if (!reverseTermExists) {
      suggestedReplacements.push(edgesReversed)
    }

    return !reverseTermExists
  })

  var fourTerms = qp.potential.filter(function(term) {
    var edges = term[1].split(",").map(e => parseInt(e))
    var edges_nonloops = edges.filter(edg => qp.edges[edg][0] !== qp.edges[edg][1])
    return edges.length === 4 && edges_nonloops.length === 4
  })

  fourTerms = fourTerms.filter(function(term) {
    //try to replace any two edges with a shortcut
    var edges = term[1].split(",").map(e => parseInt(e))
    for (var i = 0; i < edges.length; i++) {
      var thisEdge = edges[i];
      var nextEdge = edges.find(edg => qp.edges[edg][0] === qp.edges[thisEdge][1] && qp.edges[edg][1] !== qp.edges[thisEdge][0])
      if (nextEdge) {
        //can maybe be substituted with a shortcut
        var shortcut = qp.edges.findIndex(edg => edg[0] === qp.edges[thisEdge][0] && edg[1] === qp.edges[nextEdge][1])
        if (shortcut !== -1) {
          var potentialReplacement = edges.filter(edg => edg !== thisEdge && edg !== nextEdge).concat([shortcut])
          var replacementAlreadyExists = qp.potential.some(function(term) {
            return arrayEquals(term[1].split(",").map(e => parseInt(e)).sort(), potentialReplacement.sort())
          })
          
          //this term could be substituted with a 3-term
          if (!replacementAlreadyExists) {
            suggestedReplacements.push(potentialReplacement)
            return true
          }
        }
      }
    }
    return false
  })

  return {replaceable: threeTerms.concat(fourTerms), replacments: suggestedReplacements}
}
