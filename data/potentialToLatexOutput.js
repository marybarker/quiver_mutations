function remapToEdges(potential, edges) {
    return potential.map(function(term) {
        return [
            term[0], //coef,
            term[1].split(",").map(function(edgeId, idxInTerm) {
                var e = edges.find(edge => edge.id === edgeId)
                return "x_{" + e.from + "," + e.to + "}"
                /* if (idxInTerm === 0) {
                    return e.from + "," + e.to + ","
                } else if (idxInTerm === term[1].split(",").length - 1) {
                    return e.to
                } else {
                    return e.to + ","
                } */
            }).join("")
        ]
    })
}

function completeStringifyPotential(potential, edges) {
    return remapToEdges(potential, edges).map(term => term.join("")).join(" + ")
}

//eg
//4,8,6
//14,0  0,1  1,14
//14,0,1,14